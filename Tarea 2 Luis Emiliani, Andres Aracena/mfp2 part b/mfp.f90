program mfp_sampling
    !
    ! Program that samples the free flight, i.e. distance to the next interaction. and estimates its 
    ! mean value and variance for different number of histories.
    ! This program has part b (nperbatch = number of histories/trials = 0) 
    !
    ! Authors: Andres Aracena, Luis Emiliani
    ! Route: gfortran mfp.f90 mod_rng.f90 mod_mfp.f90 -o a.exe

    use mod_rng
    use mod_mfp

    implicit none

    ! Transport parameters
    integer, parameter :: ncase = 8        ! number of cases that will be studied
    integer, dimension(ncase) :: nhist     ! number of histories per case
    integer, parameter, dimension(ncase) :: nbatch = [10, 20, 30, 40, 50, 100, 200, 400] ! number of statistical batches
    integer, parameter :: nperbatch = 10   !number of histories/trials
    real, dimension(ncase,400) :: step     ! mean of each trial for each case
        ! We are using 400 like the second argument of the matrix step, because there are differents values of nbatches
        ! per case and in this problem the maximun value for nbatch is 400.
    real, dimension(ncase) :: p = 0.0      ! mean of all trials
    real, dimension(ncase) :: v = 0.0      ! variance
    real, dimension(ncase) :: r            ! relative uncertainty r = v/p

    ! Geometry parameters
    real :: sigma = 2.0                    ! total interaction cross section (cm-1)
    integer :: icase, ibatch, isample      ! loop counters
    real :: path, score                    ! auxiliary variables used for scoring

    ! Initialize the PRNG
    call rng_init(20180815)

    ! Calculation of number of histories per case
    do icase = 1, ncase
        nhist(icase) = nperbatch*nbatch(icase)
    enddo

    ! Start the sampling process. We calculate the mean and variance of the sample for 
    ! each nperbatch value given by the user.
    ncase_loop: do icase = 1, ncase
        ! Proceed with the sampling process. Start the batch loop
        batch_loop: do ibatch = 1,nbatch(icase)   
        
            ! Initialize scoring variable
            score = 0.0
            sampling_loop: do isample = 1,nperbatch 
                ! Sampling process. Accumulate the sample.
                path = mfp(sigma)
                score = score + path 
            enddo sampling_loop

            ! Accumulate results and proceed to next batch.
            step(icase,ibatch) = score/nperbatch     !mean of each trial (for each batch)

        enddo batch_loop
        
        ! Statistical analysis.

        ! Summation of the mean of each trial
        promedio: do ibatch = 1, nbatch(icase)
        p(icase) = step(icase, ibatch) + p(icase)
        enddo promedio

        ! Calculation of mean of all trials for each case
        p(icase) = p(icase)/nbatch(icase)

        ! Calculation of summation of the diference of mean of each trial and mean of all trials squared
        varianza: do ibatch = 1, nbatch(icase)
        v(icase) = (step(icase, ibatch) - p(icase))**2 + v(icase) 
        enddo varianza

        ! Calculation of sample variance s**2
        v(icase) = v(icase)/(nbatch(icase)-1)

        ! Calculation of variance
        v(icase) = sqrt(v(icase)/nbatch(icase))

        !Calculation of relative uncertainty (%)
        r(icase) = 100.0*v(icase)/p(icase)

        ! Printing of results
        write(unit=6, fmt='(I2, I5, I7, F10.5, F10.5, F10.5)') icase, nbatch(icase), nhist(icase), p(icase), v(icase), r(icase)  
    enddo ncase_loop

    ! Save results to file
    open(unit=1, file='mfp_sampling.txt')
    do icase = 1,ncase
        write(unit=1, fmt='(I2, I5, I7, F10.5, F10.5, F10.5)') icase, nbatch(icase), nhist(icase), p(icase), v(icase), r(icase)       
    enddo
    close(1)
    
end program mfp_sampling