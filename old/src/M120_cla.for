      module cla

        implicit none
        character(len=300) :: cla_fname_input = 'AA_PFM_input.TXT'
        character(len=300) :: cla_fname_gpar = 'G_PARAMETER.DAT'
        character(len=300) :: cla_fnamebase_monitor = 'walltime'
        integer(4) :: cla_monitor_onoff = 0

      contains

        subroutine input_command_line_arguments()
        implicit none

        character(len=300) :: arg  !buffer for command line arguments
        integer(4) :: i
        character(len=300) :: opt,val

        do i=1,command_argument_count()
          call get_command_argument(i, arg)
          call my_cla_parser(opt,val,arg)
          select case (trim(opt))
            case('--i')
              write(*,*)'--i given, overwrite fname_input:',trim(val)
              cla_fname_input=trim(val)
            case('--gpar')
              write(*,*)'--gpar given, overwrite fname_gpar:',trim(val)
              cla_fname_gpar=trim(val)
            case('--perfout')
              write(*,*)'--perfout given, walltime check activated:',
     &                                                         trim(val)
              cla_fnamebase_monitor=trim(val)
              cla_monitor_onoff = 1
            case default
              write(*,*)'no matching opt found. please check below:'
              call print_opt_help()
              write(*,*)'***this time, default values are used.'
          end select
        enddo !do i=1,command_argument_count()

        write(*,*)'cla_fname_input=',trim(cla_fname_input)
        write(*,*)'cla_fname_gpar=',trim(cla_fname_gpar)
        write(*,*)'cla_fnamebase_monitor=',trim(cla_fnamebase_monitor)

        end subroutine input_command_line_arguments


        subroutine my_cla_parser(opt,val,arg)
        implicit none
        character(len=300),intent(out) :: opt,val
        character(len=*),intent(in) :: arg

        character(len=300) :: arg_t
        integer(4) :: eqloc

        arg_t=arg
        arg_t=trim(arg_t)
        eqloc=scan(arg_t,'=')  !identifying the option field: find where '=' is

        opt=arg_t(1:eqloc-1)
        val=arg_t(eqloc+1:)

        end subroutine my_cla_parser

      !-------------------------------------------------------------------------------
        subroutine print_opt_help()
        implicit none

        write(*,*)'Format: [opt]=[val], no blanks allowed'
        write(*,*)'Available [opt]s and descriptions:'
        write(*,*)'  --i: Input file name, string,'
        write(*,*)'        e.g.) --i=./AA_PFM_input_AlSi.TXT'
        end subroutine print_opt_help
      end module cla
