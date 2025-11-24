! 0000000000000000000000000000000000000000000000000000000000000
! This file is part of TREKIS-4
! available at: https://github.com/N-Medvedev/TREKIS-4
! 1111111111111111111111111111111111111111111111111111111111111
! This module is written by N. Medvedev
! in 2016-2020
! 1111111111111111111111111111111111111111111111111111111111111
module Gnuplotting

implicit none 

 contains
 

subroutine process_user_gnu_parameters(gnu_extension, gnu_terminal, do_gnuplot)
   character(*), intent(in) :: gnu_extension
   character(*), intent(inout) :: gnu_terminal
   logical, intent(inout) :: do_gnuplot
   ! If user does not want gnuplot:
   if (trim(adjustl(gnu_extension)) == '0') then ! gnuplot is sad but accepts its fate
      do_gnuplot = .false.
   else ! gnuplot is happy to serve
      do_gnuplot = .true.
   endif
   ! Set terminal accordingly to the type of plots user wants:
   select case (gnu_extension)
   case ('EPS', 'Eps', 'EPs', 'eps')
      gnu_terminal = 'postscript enhanced color'
   case ('PNG', 'PNg', 'Png', 'png')
      gnu_terminal = 'pngcairo font "arial" '
   case ('JPEG', 'Jpeg', 'jpeg', 'JPG', 'Jpg', 'jpg')
      gnu_terminal = 'jpeg large font "arial"'
   case ('GIF', 'GIf', 'Gif', 'gif')
      gnu_terminal = 'gif large font "arial"'
   case ('PDF', 'PDf', 'Pdf', 'pdf')
      gnu_terminal = 'pdf color'
   end select
end subroutine process_user_gnu_parameters
 
 
subroutine write_gnu_printout(FN, first_line, last_line, file_name, x_start, x_end, y_start, y_end, col_x, col_y, lw, title, additional_info, linux_s)
   integer, intent(in) :: FN ! file nuumber to write to
   logical, intent(in) :: first_line, last_line ! is it the first line to plot? last line? for multiple curves on the same plot
   character(*) file_name ! file to plot from
   real(8), intent(in), optional :: x_start, x_end, y_start, y_end ! starting and ending points to put on axes
   character(*), intent(in), optional :: col_x, col_y ! which columns from the file to plot
   integer, intent(in), optional :: lw ! line width
   character(*), intent(in), optional :: title ! title to print on the plot
   character(*), intent(in), optional :: additional_info ! any additional gnuplot operators to pass
   logical, intent(in), optional :: linux_s ! is it linux system?
   character(500) :: format_var
   character(50) :: temp
   
   format_var = ''
   if (first_line) then 
      format_var = 'p [' 
      if (present(x_start)) then
         write(temp,'(es)') x_start
         format_var = trim(adjustl(format_var))//trim(adjustl(temp))//':'
      endif
      if (present(x_end)) then
         write(temp,'(es)') x_end
         if (present(x_start)) then
            format_var = trim(adjustl(format_var))//trim(adjustl(temp))
         else
            format_var = trim(adjustl(format_var))//':'//trim(adjustl(temp))
         endif
      endif
      format_var = trim(adjustl(format_var))//']['
      if (present(y_start)) then
         write(temp,'(es)') y_start
         format_var = trim(adjustl(format_var))//trim(adjustl(temp))//':'
      endif
      if (present(y_end)) then
         write(temp,'(es)') y_end
         if (present(y_start)) then
            format_var = trim(adjustl(format_var))//trim(adjustl(temp))
         else
            format_var = trim(adjustl(format_var))//':'//trim(adjustl(temp))
         endif
      endif
      format_var = trim(adjustl(format_var))//']'
   endif
   ! File name to get data from:
   if (present(linux_s)) then
      format_var = trim(adjustl(format_var))//' \"'//trim(adjustl(file_name))//'\" '
   else ! windows
      format_var = trim(adjustl(format_var))//' "'//trim(adjustl(file_name))//'" '
   endif
   ! Column to use for x
   format_var = trim(adjustl(format_var))//' u '//trim(adjustl(col_x))
   ! Column to use for y
   format_var = trim(adjustl(format_var))//':'//trim(adjustl(col_y))

   if (present(lw)) then ! line width
      write(temp,'(i)') lw
   else ! default value:
      write(temp,'(i)') 3
   endif
   format_var = trim(adjustl(format_var))//' w l lw '//trim(adjustl(temp))
      
   ! Any additional settings
   if (present(additional_info)) format_var = trim(adjustl(format_var))//' '//trim(adjustl(additional_info))
   
   if (present(title)) then 
      if (present(linux_s)) then
         format_var = trim(adjustl(format_var))//' title \"'//trim(adjustl(title))//'\"'
      else ! windows
         format_var = trim(adjustl(format_var))//' title "'//trim(adjustl(title))//'"'
      endif
   endif
   if (.not.last_line) format_var = trim(adjustl(format_var))//' ,\ '
   
   write(FN, '(a)') trim(adjustl(format_var))
end subroutine write_gnu_printout
 
 

subroutine write_gnuplot_script_header_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, path_sep, setkey, logx, logy)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file, gnu_term
   integer, intent(in), optional :: setkey
   logical, intent(in), optional :: logx, logy  ! set log x, or log y, scale
      
   if (present(setkey)) then
      if (path_sep .EQ. '\') then	! if it is Windows
         if ( (present(logx)) .and. (present(logy)) ) then
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey, logx, logy)
         elseif (present(logx)) then
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey=setkey, logx=logx)
         elseif (present(logy)) then
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey=setkey, logy=logy)
         else
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey=setkey)
         endif
      else ! it is linux
         if ( (present(logx)) .and. (present(logy)) ) then
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey, logx, logy)
         elseif (present(logx)) then
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey=setkey, logx=logx)
         elseif (present(logy)) then
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey=setkey, logy=logy)
         else
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey=setkey)
         endif
      endif
   else
      if (path_sep .EQ. '\') then	! if it is Windows
         if ( (present(logx)) .and. (present(logy)) ) then
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, logx=logx, logy=logy)
         elseif (present(logx)) then
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, logx=logx)
         elseif (present(logy)) then
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, logy=logy)
         else
            call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term)
         endif
      else ! it is linux
         if ( (present(logx)) .and. (present(logy)) ) then
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, logx=logx, logy=logy)
         elseif (present(logx)) then
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, logx=logx)
         elseif (present(logy)) then
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, logy=logy)
         else
            call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term)
         endif
      endif
   endif
end subroutine write_gnuplot_script_header_new


subroutine write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey, logx, logy)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(*), intent(in) :: gnu_term
   logical, intent(in), optional :: logx, logy  ! set log x, or log y, scale
   character(20) :: temp
   select case (ind)
   case(:1)	! eps
      write(FN, '(a)') '#!/bin/bash'
      write(FN, '(a)') ''
      write(FN, '(a)') 'NAME='//trim(adjustl(Out_file))
   end select
   write(FN, '(a,f3.1)') 'LW=', LW
   write(FN, '(a)') 'LABL="'//trim(adjustl(labl))//'"'
   write(temp, '(f16.3)') x_tics
!    write(FN, '(a,f6.2)') 'TICSIZ=', x_tics
   write(FN, '(a,a)') 'TICSIZ=', trim(adjustl(temp))
   write(FN, '(a)') 'echo " '
   select case (ind)
      case (:1)
!       write(FN, '(a)') 'set terminal postscript enhanced \"Helvetica\" 16 color '
      write(FN, '(a)') 'set terminal '//trim(adjustl( gnu_term ))
      write(FN, '(a)') 'set output \"$NAME\"'
      case (2:)
      write(FN, '(a)') 'set terminal x11 persist'
      write(FN, '(a)') 'unset label'
   endselect
!    write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//' \"        font \"Helvetica,20\" '
!    write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//' \"      font \"Helvetica,20\" '
   write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//' \" '
   write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//' \" '
   if (present(logx)) then
      if (logx) write(FN, '(a)') 'set logscale x'
   endif
   if (present(logy)) then
      if (logy) write(FN, '(a)') 'set logscale y'
   endif
   !write(FN, '(a)') 'set label \"$LABL\" at 150,-8 font \"Helvetica,22\" '
   if (present(setkey)) then
      select case(setkey)
      case (1)
         write(FN, '(a)') 'set key right bottom '
      case (2)
         write(FN, '(a)') 'set key left top '
      case (3)
         write(FN, '(a)') 'set key left bottom '
      case (4)
         write(FN, '(a)') 'unset key '
      case default
         write(FN, '(a)') 'set key right top '
      endselect
   else
      write(FN, '(a)') 'set key right top '
   endif
   write(FN, '(a)') 'set xtics \"$TICSIZ\" '
end subroutine write_gnuplot_script_header_linux_new



subroutine write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, gnu_term, setkey, logx, logy)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file, gnu_term
   logical, intent(in), optional :: logx, logy  ! set log x, or log y, scale
   character(20) :: temp
   select case (ind)
   case(:1)	! eps
      write(FN, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   end select
   write(FN, '(a,f3.1)') 'LW=', LW

   select case (ind)
      case (:1)
!       write(FN, '(a)') 'set terminal postscript enhanced "Helvetica" 16 color '
      write(FN, '(a)') 'set terminal '//trim(adjustl( gnu_term ))
      write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (2:)
      write(FN, '(a)') 'set terminal x11 persist'
      write(FN, '(a)') 'unset label'
   endselect
!    write(FN, '(a)') 'set xlabel "'//trim(adjustl(xlabl))//' "        font "Helvetica,20" '
!    write(FN, '(a)') 'set ylabel "'//trim(adjustl(ylabl))//' "      font "Helvetica,20" '
   write(FN, '(a)') 'set xlabel "'//trim(adjustl(xlabl))//' " '
   write(FN, '(a)') 'set ylabel "'//trim(adjustl(ylabl))//' " '
   if (present(logx)) then
      if (logx) write(FN, '(a)') 'set logscale x'
   endif
   if (present(logy)) then
      if (logy) write(FN, '(a)') 'set logscale y'
   endif
   !write(FN, '(a)') 'set label \"$LABL\" at 150,-8 font \"Helvetica,22\" '
   if (present(setkey)) then
      select case(setkey)
      case (1)
         write(FN, '(a)') 'set key right bottom '
      case (2)
         write(FN, '(a)') 'set key left top '
      case (3)
         write(FN, '(a)') 'set key left bottom '
      case (4)
         write(FN, '(a)') 'unset key '
      case default
         write(FN, '(a)') 'set key right top '
      endselect
   else
      write(FN, '(a)') 'set key right top '
   endif
   write(temp, '(f16.3)') x_tics
   write(FN, '(a,a)') 'set xtics ', trim(adjustl(temp))
end subroutine write_gnuplot_script_header_windows_new


subroutine write_gnuplot_script_ending_new(FN, File_name, path_sep)
   integer, intent(in) :: FN
   character(*), intent(in) :: File_name
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   
   if (path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      write(FN, '(a)') 'reset'
      write(FN, '(a)') '" | gnuplot '
      call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
   endif
end subroutine write_gnuplot_script_ending_new


pure subroutine cmd_vs_sh(path_sep, call_slash, sh_cmd)
   character(*), intent(in) :: path_sep
   character(*), intent(out) :: call_slash, sh_cmd
   if (path_sep .EQ. '\') then	! if it is Windows
      call_slash = 'call '
      sh_cmd = '.cmd'
   else ! It is linux
      call_slash = './'
      sh_cmd = '.sh'
   endif
end subroutine cmd_vs_sh 
 

end module Gnuplotting
