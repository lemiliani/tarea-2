program mfp_sampling
    !
    ! Program that samples the free flight, i.e. distance to the next interaction. and estimates its 
    ! mean value and variance for different number of histories.
    ! This program have parts c (nbatch = 10) and d (nbatch = 20), for diferents cases
    !
    ! Authors: Andres Aracena, Luis Emiliani
    ! Route: gfortran mfp.f90 mod_rng.f90 mod_mfp.f90 -o a.exe

    use mod_rng
    use mod_mfp

    implicit none

    ! Transport parameters
    integer, parameter :: nbatch = 20      ! number of statistical batches
    integer, parameter :: ncase = 9        ! number of cases that will be studied
    integer, dimension(ncase) :: nhist     ! number of histories per case
    integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500] !number of histories/trials
    real, dimension(ncase,nbatch) :: step  ! mean of each trial for each case
    real, dimension(ncase) :: p = 0.0      ! mean of all trials
    real, dimension(ncase) :: v = 0.0      ! variance
    real, dimension(ncase) :: r            ! relative uncertainty r = v/p

    ! Geometry parameters
    real :: sigma = 2.0                   ! total interaction cross section (cm-1)
    integer :: icase, ibatch, isample     ! loop counters
    real :: path, score                   ! auxiliary variables used for scoring

    ! Initialize the PRNG
    call rng_init(20180815)

    ! Calculation of number of histories per case
    do icase = 1, ncase
        nhist(icase) = nbatch*nperbatch(icase)
    enddo

    ! Start the sampling process. We calculate the mean and variance of the sample for 
    ! each nperbatch value given by the user.
       ncase_loop: do icase = 1, ncase
        ! Proceed with the sampling process. Start the batch loop
        batch_loop: do ibatch = 1,nbatch   
        
            ! Initialize scoring variable
            score = 0.0
            sampling_loop: do isample = 1,nperbatch(icase) 
                ! Sampling process. Accumulate the sample.
                path = mfp(sigma)
                score = score + path 
            enddo sampling_loop

            ! Accumulate results and proceed to next batch.
            step(icase,ibatch) = score/nperbatch(icase)     !mean of each trial (for each batch)

        enddo batch_loop
        
        ! Statistical analysis.

        ! Summation of the mean of each trial
        promedio: do ibatch = 1, nbatch
        p(icase) = step(icase, ibatch) + p(icase)
        enddo promedio

        ! Calculation of mean of all trials for each case
        p(icase) = p(icase)/nbatch

        ! Calculation of summation of the diference of mean of each trial and mean of all trials squared
        varianza: do ibatch = 1, nbatch
        v(icase) = (step(icase, ibatch) - p(icase))**2 + v(icase) 
        enddo varianza

        ! Calculation of sample variance s**2
        v(icase) = v(icase)/(nbatch-1)

        ! Calculation of variance
        v(icase) = sqrt(v(icase)/nbatch)

        !Calculation of relative uncertainty (%)
        r(icase) = 100.0*v(icase)/p(icase)

        ! Printing of results
        write(unit=6, fmt='(I2, I5, I7, F10.5, F10.5, F10.5)') icase, nperbatch(icase), nhist(icase), p(icase), v(icase), r(icase)  
    enddo ncase_loop    

    ! Save results to file
    open(unit=1, file='mfp_sampling.txt')
    do icase = 1,ncase
        write(unit=1, fmt='(I2, I5, I7, F10.5, F10.5, F10.5)') icase, nperbatch(icase), nhist(icase), p(icase), v(icase), r(icase)      
    enddo
    close(1)
    
end program mfp_sampling