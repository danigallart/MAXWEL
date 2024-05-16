subroutine read_solution()
    
    implicit none
    
    character*(300) :: sol_file = 'results.dat'
    integer, parameter :: sol_read_unit = 21
    character*(300) :: conn_file = 'conn_electro.dat'
    integer, parameter :: conn_read_unit = 22
    integer :: nodes,elem,i
    double precision,allocatable :: solution(:),x(:),y(:)
    integer,allocatable :: conn_electro(:,:)
    
    open(unit=sol_read_unit,file=sol_file,status='unknown')
    open(unit=conn_read_unit,file=conn_file,status='unknown')

    
    nodes = 441
    elem = 400
    
    allocate(solution(nodes),x(nodes),y(nodes))
    allocate(conn_electro(elem,4))
    
    read(sol_read_unit,*)

    
    do i=1,nodes
        read(sol_read_unit,*) x(i),y(i),solution(i)
        
    enddo
    
    
    close(sol_read_unit)
    
    read(conn_read_unit,*)

    
    do i=1,elem
        read(conn_read_unit,*) conn_electro(i,1),conn_electro(i,2),conn_electro(i,3),conn_electro(i,4)
        
    enddo
    
    
    close(conn_read_unit)



    
    
    end subroutine read_solution