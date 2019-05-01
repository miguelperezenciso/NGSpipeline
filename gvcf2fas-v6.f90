! this v should allow for haploid (two concatenated fasta) outputs
! 2018/08/06 13:28:37 : should allow for SNPs and homoz blocks colliding
! 2017-07-24 12:08:55 : allow for longer chr names
! 2016-06-21 16:08:28 : considers tetra-allelic markers (assuming several indivs)
!                       update wc function
! 2015-12-31 15:40:58 : verifies no overlapping snps and homoz blocks, discard with -nocheck option
! 2015-12-30 15:40:58 : verifies no negative blocks, discard with -nocheck option
! change for triallelic positions, allows more flexibility
! v2. rewrite to input instead of deleting reference
! sometimes even three alleles
! error when './.' genotype
! reads gvcf file from STDIN and prints to STDOUT
! WARNING: ref fa file must be in lines <= 100 characters width
!   time zcat OUT.gvcf.gz | grep -v "#" | grep -v INDEL | grep -v "\./\." | sed s/"END="// | \
!   awk -v FS="\t" -v OFS="\t" '{print $1,$2,$4,$5,$10,$8}' | cut -d";" -f1 | \
!   awk -v FS="\t" -v OFS="\t" '
!      {if( $5 ~ "0/0" )
!         print $1,$2-1,$6,$3,$3,$5
!      else
!       	 print $1,$2-1,$2,$3,$4,$5
!      }'  | sed s/':'/' '/ | sed s/"\/"/' '/ | sed s/"|"/' '/ | sed s/,//g | \
!   awk '($4!=$5 && ($6+$7)>0) || ($6==0 && $7==0)' | \
!   gvcf2fas -f ASSEMBLY.fa -nocheck > OUT.fa
!----------------
 program gvcf2fas
!----------------
 implicit none

!--> aux variables
integer          :: nchr, i, iofa=1, ios, n, g1, g2
integer(kind=8)  :: ipos, iline, chr_start(1000), nlines, line1, line2, ipos1, ipos2, last_line, &
                    last_pos, istart(1), sumbp=0, kkk
character        :: R*2, A*3, cmd*200, xc(10)*100, fafile*100=' ', chr_name(size(chr_start))*40, chr*40, &
                    chr_only=' ', last_chr*40, line*120, &
                    usage*100='USAGE: cat ID.gvcf | awk | ./gvcf2fas -f REFERENCE.fa  [-nocheck] [-hap] > ID.fa'
logical          :: check=.true., hap=.false.
integer, allocatable   :: falength(:)
character, allocatable :: faline(:)*120, gvline(:)*120, gvline2(:)*120
chr_name=' '
chr_start=0

!--> read parameters
call get_command(cmd)
call nums2(cmd, n=n, xc=xc)
do i=2, n
   select case (xc(i))
   case ('-f') !--> ref file
      fafile=xc(i+1)
   case ('-chr') !--> only this chr
      !chr_only=xc(i+1)
   case ('-hap') !--> prints two fasta files, one per haplotype
      hap=.true.
   case ('-nocheck') !--> only this chr
      check=.false.
   case ('-h')
      print*, usage; STOP
   end select
enddo
if(fafile==' ') then
   print*, USAGE; STOP
endif

!--> loads fasta in memory
nlines=wc(fafile)
allocate(faline(nlines), falength(nlines), gvline(nlines))
if(hap) allocate(gvline2(nlines))
open(iofa,file=fafile)
nchr = 0
do iline=1, nlines
   read(iofa,*) faline(iline)
   falength(iline)=len_trim(faline(iline))
   if(faline(iline)(1:1)=='>') then
      nchr=nchr+1
      if(nchr>size(chr_start)) STOP 'Increase size of chr_start'
      chr_name(nchr)=trim(faline(iline)(2:41))
      chr_start(nchr)=iline+1
   endif
enddo

!--> produces a blank reference file
line(1:maxval(falength))=repeat('N',maxval(falength))
do iline=1, nlines
   gvline(iline)=line(1:falength(iline))
enddo
if(hap) then
 do iline=1, nlines
    gvline2(iline)=gvline(iline)
 enddo
endif

!--> reads gvcf data in bed format
last_line=2
last_pos=0
last_chr = chr_name(1)
do
    read(*,*,iostat=ios) chr, ipos1, ipos2, R, A, g1, g2
    if(ios/=0) EXIT
    kkk=ipos2

    !--> ULL: assumes bed format
    ipos1=ipos1+1

    !--> NEW conditions
    !--> check if new overlapping block is homoz block, and then keep the SNP 
    if(ipos1 <= last_pos .and. chr==last_chr .and. g1==0 .and. g2==0) then
       ipos1 = last_pos+1
    endif

    !--> initializes if new chr
    if(chr/=last_chr) last_pos=0

    !--> check negative block size
    if(ipos1>ipos2 .or. ipos2 <= last_pos) then
       if(check) then
          print*, 'Negative block size!: ', chr, ipos1, ipos2
          STOP 'To ignore, use -nocheck flag'
       endif
       CYCLE
    endif
    last_chr = chr
    last_pos = ipos2

    !--> skip where no reference allele or indels, or >4 allele
    if(R=='N' .or. 'R'=='-' .or. len_trim(R)>1 .or. len_trim(A)>3) CYCLE

    !--> retrieves line and pos for each block
    istart = pack(chr_start,chr_name==trim(chr))
    line1 = istart(1) + (ipos1-1)/falength(istart(1))
    line2 = istart(1) + (ipos2-1)/falength(istart(1))
    ipos1 = ipos1 - (line1-istart(1))*falength(istart(1))
    ipos2 = ipos2 - (line2-istart(1))*falength(istart(1))

    !--> allows for triallelic SNPs: R A1,A2[,A3] must be transformed as R A1A2[A3] with sed
    if(len_trim(A)==2 .or. len_trim(A)==3) then
       if(g1==g2 .and. g1>0) then      !--> hom for non ref allele
          A(1:1)=A(g1:g1); g1=1; g2=1
       elseif(g1/=g2) then           
          if(g1>0 .and. g2>0) then     !--> het for both non ref alleles
             R(1:1)=A(g1:g1); A(1:1)=A(g2:g2); g1=0; g2=1
          elseif(g1==0) then  !--> g1=R, to output phased haploids
             !g2=max(g1,g2)
             A(1:1)=A(g2:g2); g1=0; g2=1
          elseif(g2==0) then  !--> g2=R
             !g2=max(g1,g2)
             A(1:1)=A(g1:g1); g1=1; g2=0
          else
             STOP 'You should not have reached this'
          endif
       endif
    endif

    !--> fillin genotypes or reference alleles
    select case(g1 + g2)
     case(0)
        if(line1==line2) then
           do i=ipos1,ipos2; gvline(line1)(i:i)=faline(line1)(i:i); enddo
        else
           do i=ipos1, falength(line1); gvline(line1)(i:i)=faline(line1)(i:i); enddo
           do iline=line1+1, line2-1
              do i=1, falength(iline); gvline(iline)(i:i)=faline(iline)(i:i); enddo
           enddo
           do i=1, ipos2; gvline(line2)(i:i)=faline(line2)(i:i); enddo
        endif
     case(1)
       gvline(line2)(ipos2:ipos2)=get_iupac(R(1:1)//A(1:1))
     case(2)
       gvline(line2)(ipos2:ipos2)=A(1:1)
    end select

    !--> if hap, only needs to change het positions
    if(hap) then
    select case(10*g1 + g2)
     case(0)
        if(line1==line2) then
           do i=ipos1,ipos2; gvline2(line1)(i:i)=faline(line1)(i:i); enddo
        else
           do i=ipos1, falength(line1); gvline2(line1)(i:i)=faline(line1)(i:i); enddo
           do iline=line1+1, line2-1
              do i=1, falength(iline); gvline2(iline)(i:i)=faline(iline)(i:i); enddo
           enddo
           do i=1, ipos2; gvline2(line2)(i:i)=faline(line2)(i:i); enddo
        endif
     case(1)
       gvline(line2)(ipos2:ipos2) =R(1:1)
       gvline2(line2)(ipos2:ipos2)=A(1:1)
     case(10)
       gvline(line2)(ipos2:ipos2) =A(1:1)
       gvline2(line2)(ipos2:ipos2)=R(1:1)
     case(11)
       gvline2(line2)(ipos2:ipos2)=A(1:1)
    end select
    endif
!print*, kkk, g1, g2, R, A, gvline(line2)(ipos2:ipos2)
enddo

!--> rewrite chr names just in case have been overwritten
if(hap) then
   do i=1, nchr
      gvline(chr_start(i)-1) ='>'//trim(chr_name(i))//' h1'
      gvline2(chr_start(i)-1)='>'//trim(chr_name(i))//' h2'
   enddo
else
   do i=1, nchr
      gvline(chr_start(i)-1)='>'//trim(chr_name(i))
   enddo
endif

do iline=1, nlines
   write(*,'(a)') trim(gvline(iline))
enddo
if(hap) then
 do iline=1, nlines
   write(*,'(a)') trim(gvline2(iline))
 enddo
endif

CONTAINS

!------------------
 function get_iupac (bp)
!------------------
 character :: get_iupac*1, bp*(*), snp(17)*2, c(1)*1, iupac(17)*1
 iupac = (/'A',  'C',  'T',  'G',  'R',  'R',  'Y',  'Y',  'K',  'K',  'M', 'M',  'S',  'S',  'W', 'W',  'N'/)
 snp =   (/'AA', 'CC', 'TT', 'GG', 'AG', 'GA', 'TC', 'CT', 'GT', 'TG', 'AC','CA', 'CG', 'GC', 'AT','TA', 'NN'/)
 c = pack(iupac, snp==trim(bp))
 get_iupac = c(1)
!------------
 end function
!------------

!----------------
 INTEGER FUNCTION wc(infile)
!----------------
! counts file lines; removes blank end line if it exists
 character :: infile*(*), ofile*20
 integer   :: ioin=333, ios
 integer(8):: time(8)
 logical   :: ex
 call date_and_time(values=time)
 time(1)=10000*time(8)+sum(time)
 write(ofile,'(i12)') time(1)
 ofile=trim(adjustl(ofile))//'.tzx'

 wc = 0
 inquire(file=infile,exist=ex)
 !--> all this fuss is to avoid problems in very wide files, rm trailing blanks and rm last empty line !!
 if(ex) then
    call system("awk '{print substr($1,1,10)}' " // trim(adjustl(infile)) //  &
                " | awk 'NR > 1{print t} {t = $0}END{if (NF) print }' | wc -l > "//ofile)
    open(ioin,file=ofile)
    read(ioin,*,iostat=ios) wc
    if(ios/=0) wc=0
    close (ioin, status='delete')
 endif
!------------
 END FUNCTION
!------------

!----------------
 subroutine nums2 (a, n, x, xc)
!--------------------------------------------------------------------------
! separates array a into items delimited by blanks. character elements are
! put into optional character vector xc, decoded numeric values
! into optional real vector x, and n contains the number of items. The
! dimension of x and xc can be lower than n.
! format changed slightly to read real numbers in scientific notation (mpe)
! 2/23/2000 IMisztal
!--------------------------------------------------------------------------

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
         if (numeric(adjustl(a(first:last)))) then    !NEW (mpe)
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
! error corrected in init z
 character :: a*(*)
 integer   :: i, z
 numeric = .false.
 z=0
 do i=1, len_trim(a)
    if(scan(a(i:i),'.1234567890')==0) RETURN
    if(a(i:i)=='.') z=z+1
 enddo
 if(z<=1) numeric = .true.
!------------
 end function
!------------

!-----------
 END PROGRAM
!-----------
