      subroutine opencount(maxev)
      implicit none
      include 'pwhg_rnd.h'
      integer maxev
      character * 30 file
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ios
      character * 7 string
      real * 8 powheginput
      external powheginput
      integer nev,j,iun,iret
      common/copencount/iun
      maxev=0
      file='pwgevents.lhe'
      call pwhg_io_open_read(trim(file),iun,ios)
c      open(unit=iun,file=file,status='old',iostat=ios)
      if(ios.ne.0) then
         do j=30,0,-1
            if(file(j:j).ne.' ') exit
         enddo
         write(*,*)' file not found:',file(1:j)
         write(*,*)' enter name of event file'
         read(*,'(a)') file
         call pwhg_io_open_read(trim(file),iun,ios)
c         open(unit=iun,file=file,status='old',iostat=ios)
         if(ios.ne.0) then
            write(*,*) 'cannot open; aborting ...'
            call exit(-1)
         endif
c get the name prefix
         j=index(file,'events')
         if(j.gt.1) then
            pwgprefix=file(1:j-1)
            lprefix=j-1
         else
            lprefix=0
         endif
c or a sequence number for manyseeds usage
         if(file(j+6:j+6).eq.'-') then
            rnd_cwhichseed=file(j+7:j+10)
         else
            rnd_cwhichseed='none'
         endif
      else
         pwgprefix='pwg'
         lprefix=3
         rnd_cwhichseed='none'
      endif
      write(*,*) 'prefix: "'//pwgprefix(1:lprefix)//'",', ' seed: "'//
     1     rnd_cwhichseed//'"'
c
      write(*,*) ' Opened event file ',file
      write(*,*) ' Counting events in ', file
      write(*,*) ' This may take some time...'
 1    continue
      call pwhg_io_read(iun,string,iret)
      if(iret /= 0) goto 2
c      read(unit=iun,fmt='(a)',end=2) string
      if(string.eq.'</event') then
         maxev=maxev+1
         goto 1
      endif
      goto 1
 2    continue
      write(*,*) ' Found ',maxev,' events in file ',file
      if (maxev.eq.0) then
         write(*,*) ' NO EVENTS!! Program exits'
         call exit(3)
      endif
      call pwhg_io_rewind(iun)
      end

      subroutine opencountunit(maxev,iun)
      implicit none
      include 'pwhg_rnd.h'
      integer maxev,iun
      character * 30 file
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer ios
      character * 7 string
      real * 8 powheginput
      external powheginput
      integer nev,j
      call newunit(iun)
      maxev=0
      file='pwgevents.lhe'
      call pwhg_io_open_read(trim(file),iun,ios)
c      open(unit=iun,file=file,status='old',iostat=ios)
      if(ios.ne.0) then
         do j=30,0,-1
            if(file(j:j).ne.' ') exit
         enddo
         write(*,*)' file not found:',file(1:j)
         write(*,*)' enter name of event file'
         read(*,'(a)') file
         call pwhg_io_open_read(trim(file),iun,ios)
c         open(unit=iun,file=file,status='old',iostat=ios)
         if(ios.ne.0) then
            write(*,*) 'cannot open; aborting ...'
            call exit(-1)
         endif
c get the name prefix
         j=index(file,'events')
         if(j.gt.1) then
            pwgprefix=file(1:j-1)
            lprefix=j-1
         else
            lprefix=0
         endif
c or a sequence number for manyseeds usage
         if(file(j+6:j+6).eq.'-') then
            rnd_cwhichseed=file(j+7:j+10)
         else
            rnd_cwhichseed='none'
         endif
      else
         pwgprefix='pwg'
         lprefix=3
         rnd_cwhichseed='none'
      endif
      write(*,*) 'prefix: "'//pwgprefix(1:lprefix)//'",', ' seed: "'//
     1     rnd_cwhichseed//'"'
c
      write(*,*) ' Opened event file ',file
      write(*,*) ' Counting events in ', file
      write(*,*) ' This may take some time...'
 1    continue
      call pwhg_io_read(iun,string,ios)
      if(ios /= 0) goto 2
c     read(unit=iun,fmt='(a)',end=2) string
      if(string.eq.'</event') then
         maxev=maxev+1
         goto 1
      endif
      goto 1
 2    continue
      write(*,*) ' Found ',maxev,' events in file ',file
      if (maxev.eq.0) then
         write(*,*) ' NO EVENTS!! Program exits'
         call exit(3)
      endif
      call pwhg_io_rewind(iun)
      end
