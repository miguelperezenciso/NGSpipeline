! small change in header with =
! prints n missing optionally
! output in format tped012
! program to convert fasta in vcf simplified format for a set of aligned files
! all should have the same width in fasta file
! vcf2fas -i fasta.list  -f reference.fa | sed s/' '//g > out.vcf
! MPE: miguel.perez@uab.es
! #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002s

!----------------
 program fas2vcf
!----------------
!use aux_sub

implicit none
integer                :: ioref=10, i, io=1, iout=7, iofam=8, iomiss=2, ios, n, nfiles, ibp, winsize=0, iwin, &
                          current_iwin
integer, allocatable   :: iofa(:), miss(:)
integer(kind=8)        :: ipos, bplength
real                   :: x(20)
logical                :: ex, tped=.false., c012=.false., printmiss=.false.
integer(kind=8), allocatable :: lmiss(:)
character, allocatable       :: line(:)*100, bp(:)*1, vcf(:)*3, name(:)*50
character                    :: falist*100='', faref*200='', outfile*100='fas2vcf.out', faline*100, chr*10=' ', &
                                cmd*200, xc(size(x))*200, bpr*1, c*100, alleles*5, region*100=' ', usage*200, &
                                missfile*30='fas2vcf.miss', current_chr*20
usage = "USAGE: ./fas2vcf -i fasta.list -f reference.fa | sed s/' '//g | gzip -c > out.vcf.gz"
call get_command (cmd)
call nums2 (cmd, n, x, xc)
do i=2, n
   select case (xc(i))
   case ('-h', '-help')
      print*, trim(usage)
   case ('-i')
      falist=xc(i+1)
   case ('-f')
     faref=xc(i+1)
     inquire(file=faref,exist=ex)
     if(.not.ex) then
        print*, trim(usage)
        STOP 'Reference fasta file required'
     endif
   case ('-tped')
     tped=.true.
   case ('-012')
     tped=.true.
     c012=.true.
   case ('-printmiss')
     !missfile=xc(i+1)
     open(iomiss,file=missfile)
     winsize=100000
     printmiss=.true.
   case ('-w')
     winsize=x(i+1)
     if(winsize<1) STOP 'Window size must be positive'
   end select
enddo
if(faref=='' .or. falist=='') then
   print*, trim(usage)
   STOP 'Reference fasta and file with fasta files must be specified with -f and -i'
endif

!--> count and open fasta files
nfiles=wc(falist)
if(nfiles==0) STOP 'File with fasta files and sample names is required, one per line'
allocate(line(nfiles), iofa(nfiles), vcf(nfiles), name(nfiles), lmiss(nfiles), miss(nfiles) )
iofa=[(ioref+i, i=1, nfiles)]
open(io,file=falist)
do i=1, nfiles
   read(io,*) xc(1), name(i)
   open(iofa(i),file=xc(1))
enddo

!--> open reference fasta
open(ioref,file=faref)

!---> open outfiles
open(iout,file='qq')

!--> writes header
! #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002s
write(*,*) '##fileformat=VCFv4.1'
write(*,*) '##program=fas2vcf'
do i=1, nfiles
   write(*,*) '##samples=',trim(name(i))
enddo
write(*,*) '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
write(*,*) '#CHROM',achar(9),'POS',achar(9),'ID',achar(9),'REF',achar(9),'ALT',achar(9),'QUAL',&
              achar(9),'FILTER',achar(9),'INFO',achar(9),'FORMAT',(achar(9),trim(name(i)),i=1,nfiles)


!--> main loop, reads through different chromosomes and fasta files
ipos=0
lmiss=0
current_iwin=1
do
   read(ioref,'(a)',iostat=ios) faline
   if(ios/=0) EXIT

   !--> new chr
   if(faline(1:1)=='>') then
      if(ipos>0 .and. printmiss) write(iomiss,*) current_chr, current_iwin, bplength, lmiss(:)
      ipos=0
      bplength=0
      lmiss=0
      call nums2 (faline, n, x, xc)
      chr=trim(adjustl(xc(1)(2:len_trim(xc(1)))))
      do i=1, nfiles
         read(iofa(i),'(a)') xc(1)
         c=trim(adjustl(xc(1)(2:len_trim(xc(1)))))
         if(c.ne.chr) then
            print*, 'in unit ',i,' chrs not match', chr,' ',c
            STOP
         endif
      enddo
      current_chr=chr
      CYCLE !--> next line of reference if new chr
   endif

   !--> read fasta line for each sample
   do i=1, nfiles
      read(iofa(i),*) line(i)
      if(len_trim(line(i))/=len_trim(faline)) STOP 'Fasta lines must be of equal length'
   enddo

   !--> process sequence line for every sample
   do ibp=1, len_trim(adjustl(faline))
      ipos=ipos+1
      alleles=''
      !--> dirty work is left to subroutine, if reference allele is known
      if(faline(ibp:ibp)/='N') then
         call get_vcf(faline(ibp:ibp), line(:)(ibp:ibp), alleles, vcf, miss)
         bplength=bplength+1
         lmiss=lmiss+miss
      endif
      !--> A SNP exists
      if(alleles.ne.'') then
          if(.not.tped) then
             write(*,*) trim(chr),achar(9), ipos,achar(9), '.',achar(9), faline(ibp:ibp),achar(9), &
                           trim(alleles),achar(9), '30',achar(9), 'PASS',achar(9), '.',achar(9),  &
                          'GT',  (achar(9),vcf(i),i=1,size(vcf))
          elseif (tped .and. len_trim(alleles)==1) then
             write(cmd,*) ipos
             write(*,*) chr, trim(chr)//'_'//trim(adjustl(cmd)), 0, ipos, (' ',vcf2_012(vcf(i)),i=1,size(vcf))
          endif
      endif
      !--> printmiss
      if(printmiss) then
         iwin = ipos/winsize+1
         if(iwin .ne. current_iwin) then
            write(iomiss,*) current_chr,' ', current_iwin, bplength, lmiss(:)
            current_iwin = iwin
            lmiss = 0
            bplength = 0
         endif
      endif
   enddo

enddo
if(printmiss) write(iomiss,*) current_chr, iwin, bplength, lmiss(:)

CONTAINS

!------------------
 subroutine get_vcf (a,ai,alleles,vcf,miss)
!------------------
 character              :: a*(*), ai(:)*(*), alleles*(*), vcf(:)*(*),  c*1
 character, allocatable :: genotypes(:)*2
 integer                :: i, h, n, iallele, miss(:)
 alleles=' '
 n=size(ai)
 vcf=''
 miss=0

 do i=1, n
    if(ai(i)=='N' .or. ai(i)=='.') miss(i)=1
 enddo

 !--> check whether a potential SNP
 if(.not.all(a==pack(ai,ai/='N'))) then
    select case (a)
    case('A','C','G','T')

    allocate (genotypes(n))

    !--> returns genotype
    do i=1, n
       genotypes(i)=get_snp(ai(i))
    enddo
!print*, 'gp',ai,' ',genotypes

    !--> returns n alleles
    call parse_alleles(alleles,a,genotypes)
!print*, a,'alleles',alleles

    !--> genotype scanned
    do i=1, n
       do h=1,2
          if(genotypes(i)(h:h)==a) then
             vcf(i)=trim(vcf(i))//'0'
          else
             do iallele=1,len_trim(alleles)
                if(genotypes(i)(h:h)==alleles(iallele:iallele)) then
                   write(c,'(i1)') iallele
                   vcf(i)=trim(vcf(i))//c
                endif
             enddo
          endif
       enddo
       !--> insert '/' or set to missing
       if(len_trim(vcf(i))>0) then
          vcf(i) = reorder_vcf(vcf(i))
          vcf(i)(3:3)=vcf(i)(2:2)
          vcf(i)(2:2)='/'
       else
          vcf(i)='./.'
       endif
    enddo

    !--> transform alleles AT into A,T ..
    if(len_trim(alleles)==2) then
       alleles(3:3)=alleles(2:2)
       alleles(2:2)=','
    elseif(len_trim(alleles)==3) then
       alleles(5:5)=alleles(3:3)
       alleles(3:3)=alleles(2:2)
       alleles(2:2)=','
       alleles(4:4)=','
    endif
    deallocate(genotypes)

    end select
 endif
!-------------
 end subroutine
!-------------

!------------------------
 subroutine parse_alleles (alleles,a,gi)
!------------------------
! returns index of non reference alleles
 character :: a*(*), gi(:)*(*), alleles*(*)
 integer   :: i, j, freq(4)
 alleles=''
 freq=0

 do j=1, size(gi)
    do i=1,2
      select case (gi(j)(i:i))
       case('A','a')
         if(freq(1)==0.and.a/='A') alleles=trim(adjustl(alleles))//'A'
         freq(1)=freq(1)+1
       case('C','c')
         if(freq(2)==0.and.a/='C') alleles=trim(adjustl(alleles))//'C'
         freq(2)=freq(2)+1
       case('G','g')
         if(freq(3)==0.and.a/='G') alleles=trim(adjustl(alleles))//'G'
         freq(3)=freq(3)+1
       case('T','t')
         if(freq(4)==0.and.a/='T') alleles=trim(adjustl(alleles))//'T'
         freq(4)=freq(4)+1
      end select
    enddo
 enddo
!--------------
 end subroutine
!--------------

!----------------
 function get_snp (bp)
!----------------
 character :: get_snp*2, bp*(*), snp(12)*2, c(1)*2, iupac(12)*1
 iupac = (/'A',  'C',  'T',  'G',  'R',  'Y',  'K',  'M',  'S',  'W',  'N',  '-'/)
 snp =   (/'AA', 'CC', 'TT', 'GG', 'AG', 'CT', 'GT', 'AC', 'GC', 'AT', 'NN', 'NN'/)
 c = pack(snp, iupac==trim(bp))
 get_snp = c(1)
!----------------
 end function
!----------------

!--------------------
 function reorder_vcf (vcf)
!--------------------
 character :: reorder_vcf*2, vcf*(*), vcf0(16)*2, c(1)*2, vcf1(16)*2
 vcf0 = (/'00',  '01',  '02',  '03',  '10',  '11',  '12',  '13',  '20',  '21',  '22', '23', '30', '31', '32', '33'/)
 vcf1 = (/'00',  '01',  '02',  '03',  '01',  '11',  '12',  '13',  '02',  '12',  '22', '23', '03', '13', '23', '33'/)
 c = pack(vcf1, vcf0==trim(vcf))
 reorder_vcf = c(1)
!----------------
 end function
!----------------

!-------------------------
 integer function vcf2_012 (vcf)
!-------------------------
 integer   :: icode(4)=[9,0,1,2], ic(1)
 character :: vcf*(*), vcode(4)*3=['./.', '0/0', '0/1', '1/1']
 ic = pack(icode, vcode==trim(vcf))
 vcf2_012 = ic(1)
!----------------
 end function
!----------------


!-----------------------
 subroutine parse_region (region, chr, r_pos1, r_pos2)
!-----------------------
 character*(*)    :: region, chr
 integer (kind=8) :: r_pos1, r_pos2, i1, ic, l
 chr=region
 ic = index(region,':')
 i1 = index(region,'-')
 l = len_trim(region)
 if (ic > 0) then
    read(region(1:ic),*) chr
    if(i1==0) then
       if (ic+1<l) read(region(ic+1:l),*) r_pos1
    else
       if(ic+1<i1-1) read(region(ic+1:i1-1),*) r_pos1
       if(i1+1<l) read(region(i1+1:l),*) r_pos2
    endif
 endif
!--------------
 end subroutine
!--------------

!----------------
 INTEGER FUNCTION wc(infile)
!----------------
 character :: infile*(*), ofile*3='tzx'
 integer   :: ioin=333, ios
 logical   :: ex
 wc = 0
 inquire(file=infile,exist=ex)
 if(ex) then
    call system('wc -l '//infile//'>'//ofile)
    open(ioin,file=ofile)
    read(ioin,*) wc
    close (ioin, status='delete')
 endif
!------------
 END FUNCTION
!------------


!------------------
function tab2space (a)
!------------------
! replace all TAB with SPACE (I. Aguilar)

character*(*) :: a
character(len=len(a)) tab2space
integer :: i

do i=len_trim(a),1,-1
if (a(i:i)==achar(9)) then
a(i:i) = " "
endif
enddo
tab2space=a
!------------
end function
!------------

!----------------
subroutine nums2 (a, n, x, xc)
!----------------
! separates array a into items delimited by blanks or tabs.
! character elements are put into optional character vector xc, decoded numeric values
! into optional real vector x, and n contains the number of items. The
! dimension of x and xc can be lower than n.
! format changed slightly to read real numbers in scientific notation (mpe)
! 2/23/2000 IMisztal

character (*)          :: a
character (*),optional :: xc(:)
real,optional          :: x(:)
integer :: n, curr, first, last, lena, stat, i

curr=1
lena=len(a)
n=0

!do i=1, lena
!   if(ichar(a(i:i))==9) a(i:i)=" "
!enddo

do
!--> search for first nonspace
first=0
do i=curr,lena
if (a(i:i) /= ' ' ) then
first=i
exit
endif
enddo
if (first == 0) exit

!--> search for first space
curr=first+1
last=0
do i=curr,lena
if (a(i:i) == ' ' ) then
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
print*, a,n,'NUMS2: Increase size of XC'
endif
endif

if (present(x)) then
if (size(x) >= n) then
if (numeric(a(first:last))) then    !NEW (mpe)
read(a(first:last),*) x(n)
else
x(n)=0
endif
else
print*, 'NUMS2: Increase size of X'
endif
endif

curr=last+1
enddo
!--------------
end subroutine
!--------------


!------------------------
 logical function numeric (a)
!------------------------
 character :: a*(*)
 integer :: i, z=0
 numeric = .false.
 do i=1, len_trim(a)
    if(scan(a(i:i),'.1234567890')==0) RETURN
    if(a(i:i)=='.') z=z+1
 enddo
 if(z<=1) numeric = .true.
!------------
 end function
!------------


!------------------
 function get_iupac (bp)
!------------------
 character :: get_iupac*1, bp*(*), snp(26)*2, c(1)*1, iupac(26)*1
 iupac = (/'A',  'C',  'T',  'G',  'A',  'C',  'T',  'G',  'R',  'Y',  'K',  'M',  'S',  'W',  'R',  'Y',  'K',  'M',  'S',  'W',  &
                                                           'R',  'Y',   'K',   'M',   'S',   'W'/)
 snp =   (/'A ', 'C ', 'T ', 'G ', 'AA', 'CC', 'TT', 'GG', 'AG', 'CT', 'GT', 'AC', 'GC', 'AT', 'GA', 'TC', 'TG', 'CA', 'CG', 'TA', &
                                                           'R ', 'Y ',  'K ',  'M ',  'S ',  'W '/)
 c = pack(iupac, snp==bp)
 get_iupac = c(1)
!----------------
 end function
!----------------

end program
