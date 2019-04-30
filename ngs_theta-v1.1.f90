!----------------
 subroutine usage
!----------------
 print'(50(a/))', &

' LAST MODIFICATIONS',&
' Tue 22 Oct 2013 : small change to allocate phf if maxrd<15',&
' Wed 13 Feb 2013 10:14:54 CET: rewriting, piping and line commands allowed',&
' Mon 14 Jan 2013 10:13:40 CET: uses Lynch coeff for > 15 depth in phf',&
'',&
' NGS_THETA v. 1',&
'   Implements Lynch / Ferretti estimator of theta from single individual sequencing for low coverage',&
'   as described in ',&
'     Partial short-read sequencing of a highly inbred Iberian pig and genomics inference thereof',&
'     Esteve-Codina et al. 2011, Heredity, 107:256-64. doi: 10.1038/hdy.2011.13',&
'     http://www.nature.com/hdy/journal/v107/n3/full/hdy201113a.html',&
'   To be used with covXwin.pl',&
'',&
' USAGE',&
'    samtools mpileup -Bq $MAPQ -d 100000 $SAMPLE.bam  | \' ,&
'    perl covXwin.pl -v $SAMPLE.vcf -w $WINSIZE -d $MINCOV -m $MAXCOV -b $SAMPLE.bam | \ ' ,&
'    ngs_theta -d $MINCOV -m $MAXCOV > $OUT.wintheta ',&
'',&
'    $MAPQ is minimum RMS quality [20]',&
'    $SAMPLE is sample name (req)',&
'    $WINSIZE is window size in bp [100000]',&
'    $MINCOV is minimum depth [req]',&
'    $MAXCOV is maximum depth [req]',&
'',&
' NOTES',&
'    To compile: f95 fact-m.f90 ngs_theta-v1.1.f90 -o ngs_theta',&
'    The input (from covXwin.pl) format is, one line per window ',&
'       CHR.window, S, F, I, cov(min:max)',&
'    Output format is',&
'       CHR   WINDOW NSNP   NFIX   BP  DEPTH  THETA_HET  THETA_FIX  SNP_HAT  FIX_HAT  THETA_HKA  CHI2',&
'    where',&
'       NSNP:      raw number of heterozygous positions',&
'       NFIX:      raw number of fixed differences',&
'       BP:        bp covered in that window, filtered by restrictions min & amx depth and RMS quality',&
'       DEPTH:     average depth',&
'       THETA_HET: estimated heterozygosity (eq. 3 in Esteve-Codina) ',&
'       THETA_FIX: estimated divergence ',&
'       SNP_HAT:   adjusted no. of snps (eq. 4)',&
'       FIX_HAT:   adjusted no. of fixed differences (eq. 5)',&
'       THETA_HKA: HKA estimate of diversity (SNP-hat + F-hat)/(T+2) in Esteve-Codina',&
'       CHI2:      Chis2 of hka test',&
' WARNING!',&
'    Utmost care in interpreting HKA as is not outgroup!!!!'
 STOP
!--------------
 end subroutine
!--------------

!----------------
 program ngstheta
!----------------
 use factorialmodule
 
 implicit none
 integer :: minrd=0, maxrd=0, rdbar=13, lwindow=500000
 integer :: iout=10, maxwin, rd, ichr, ipos, iwin, i, j, k, &
            nsnp, ndel, io, nwin=0, ifile, n, nfix, kmax
 integer, allocatable :: winsnp(:,:), winchr(:,:), snpos(:,:), work(:), winfix(:,:), &
                         fpos(:,:), winlength(:)
 real, allocatable :: ph(:), pf(:), phf(:)
 real :: q=0, theta, cov, fheta, snp_hat, fix_hat, t, vs, vf, hka, chi2, theta_hka, &
         x(100)
 real (selected_real_kind(15, 307)) :: es=0, ef=0
 character :: filesnp*100, a*1, chr*30, xc(size(x))*100, cmd*200

 call get_command(cmd)
 call nums2(cmd, n, x, xc)
 do i=2, n
   select case (xc(i))
     case ('-d')
        minrd = x(i+1)
     case ('-m')
        maxrd = x(i+1)
     case ('-h')
        call usage
   end select
 enddo
 if (minrd == 0 .or. maxrd == 0) call usage
 
!--> now computes pH(nr) ,ie, empirical error
 allocate (ph(maxrd), phf(maxrd), pf(maxrd), winlength(maxrd))
 
 !-- ferretti's coefficient (mac > 0.2 and two copies non ref allele)
 phf=0
 do i=minrd, min(15,maxrd)
    phf(i)=1.
    kmax = max(1,(i-1)/5)
    do j=0, kmax
       phf(i) = phf(i) - combinatorial(i,j) * 0.5**(i)
    enddo
    kmax = 1
    do j=i-(kmax-1),i
       phf(i) = phf(i) - combinatorial(i,j) * 0.5**(i)
    enddo
 enddo
 !--> Lynch's coefficient for enough depth
 if(maxrd>=15) then
    do i=15, maxrd
       phf(i) = 1 - 0.5**(i-1)
    enddo
 endif

 winlength=0
 do
    read(*,*,iostat=io) chr, iwin, nsnp, nfix, ndel, winlength(minrd:maxrd)
    if(io/=0) EXIT
    nwin=nwin+1
    if(sum(winlength(:))>0) then
       theta = nsnp / dot_product(winlength(:), phf)
       snp_hat = theta * sum(winlength(:))
       pf = (/(theta*0.5**j, j=1, maxrd)/) !-> P of not being a fixed difference
       pf(1:minrd-1)=0
       fix_hat = max(0.0, nfix - dot_product(winlength(:), pf)) 
       fheta = fix_hat / sum(winlength(:))
       cov = dot_product(winlength(:), (/(j, j=1,maxrd)/) ) / real(sum(winlength(:)))
    else
       cov = 0; theta = 0; fheta=0; snp_hat=0; fix_hat=0
    endif
    es = es + snp_hat
    ef = ef + fix_hat
    write(iout,*) chr,' ', iwin, nsnp, nfix, sum(winlength(:)), cov, theta, fheta, snp_hat, fix_hat
 enddo
 if(nwin==0) call usage
 close (iout)

 !--> HKA
 open(iout)
 t = ef/es-1.
 print*, 'CHR   WIN NSNP   NFIX   BP  DEPTH  THETA_HET  THETA_FIX  SNP_HAT  FIX_HAT  THETA_HKA  CHI2'
 do i=1, nwin
    read(iout,*) chr, iwin, nsnp, nfix, k, cov, theta, fheta, snp_hat, fix_hat
    es = (snp_hat + fix_hat) / (t+2.)
    ef = es * (t+1.)
    vs = es*(es+1.)
    vf = ef*(ef+1)
    if (vs>0 .and. k>0) then
       chi2 = (snp_hat-es)**2/vs + (fix_hat-ef)**2/vf
       theta_hka=real(es/k)
    else
       chi2 = 0.
       theta_hka=0
    endif
    write(*,*) adjustl(trim(chr)),' ', iwin, nsnp, nfix, k, cov, theta, fheta, snp_hat, fix_hat, theta_hka, chi2
 enddo 
 close(iout, status='delete')

CONTAINS

!----------------
 subroutine nums2 (a, n, x, xc)
!----------------
! separates array a into items delimited by blanks. character elements are
! put into optional character vector xc, decoded numeric values 
! into optional real vector x, and n contains the number of items. The 
! dimension of x and xc can be lower than n.
! A modification of nums() from f77 to f90
! Now accepts real numbers
! format changed slightly to read real numbers in scientific notation (mpe)
! 2/23/2000 IMisztal

 character (*)          :: a
 character (*),optional :: xc(:)
 real,optional          :: x(:)
 integer :: n, curr, first, last, lena, stat, i

 curr=1
 lena=len(a)
 n=0

 do 
   !--> search for first nonspace
   first=0
   do i=curr,lena
     if (a(i:i) /= ' ') then
        first=i
        exit
     endif
   enddo
   if (first == 0) exit
  
   !--> search for first space
   curr=first+1
   last=0
   do i=curr,lena
      if (a(i:i) == ' ') then
        last=i
        exit
      endif
   enddo
  
   if (last == 0) last=lena
  
   n=n+1
   if (present(xc)) then
      if (size(xc) >= n) then
         xc(n)=a(first:last)
      else
         print*, 'NUMS2: Increase size of XC'
      endif
   endif   

   if (present(x)) then
      if (size(x) >= n) then
        read(a(first:last),*,iostat=stat) x(n)    !NEW (mpe)
        if (stat /=0) x(n)=0
      else
        print*,  'NUMS2: Increase size of X'       
      endif
   endif   
 
   curr=last+1
 enddo
!--------------
 end subroutine
!--------------

!-----------
 end program
!-----------
