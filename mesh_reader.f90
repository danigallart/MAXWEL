subroutine mesh_reader
    use def_io
    use def_variables
    use def_vectors
    
    implicit none
    
    integer :: istat, ii
    
    ! Node's coordinates
    
    READ(nodes_unit,*,IOSTAT=istat) NP
    
    allocate(coorx(NP), coory(NP))
    do ii=1,NP
        read(nodes_unit,*,iostat=istat) coorx(ii),coory(ii)
    end do
    
    CLOSE(nodes_unit)
        
    ! Conectivity matrix
    
    nodpel = 3
    
    read(elements_unit,*,iostat=istat) NE
    
    allocate(conectividad(NE,nodpel))
    DO ii=1,NE
        read(elements_unit,*,iostat=istat)  conectividad(ii,nodpel-2), conectividad(ii,nodpel-1), conectividad(ii,nodpel)
    END DO
    
    CLOSE(elements_unit)
    
    ! PML region and boundary nodes
    
    read(pml_bin_nodes_unit,*,iostat=istat) n_pml_bin
    
    allocate(pml_bin_nodes(n_pml_bin))
    
    do ii=1,n_pml_bin
        read(pml_bin_nodes_unit,*,iostat=istat) pml_bin_nodes(ii)
    end do
    
    close(pml_bin_nodes_unit)
    
    read(pml_bout_nodes_unit,*,iostat=istat) n_pml_bout
    
    allocate(pml_bout_nodes(n_pml_bout))
    
    do ii=1,n_pml_bout
        read(pml_bout_nodes_unit,*,iostat=istat) pml_bout_nodes(ii)
    end do
    
    close(pml_bout_nodes_unit)
    
    read(pml_nodes_unit,*,iostat=istat) n_pml
    
    allocate(pml_nodes(n_pml))
    
    do ii=1,n_pml
        read(pml_nodes_unit,*,iostat=istat) pml_nodes(ii)
    end do
    
    close(pml_nodes_unit)
    
    ! Scattering region and boundaries
    
    read(scatb_nodes_unit,*,iostat=istat) n_scatb
    
    allocate(scatb_nodes(n_scatb))
    
    do ii=1,n_scatb
        read(scatb_nodes_unit,*,iostat=istat) scatb_nodes(ii)
    end do
    
    close(scatb_nodes_unit)
    
    read(scatin_nodes_unit,*,iostat=istat) n_scatin
    
    allocate(scatin_nodes(n_scatin))
    
    do ii=1,n_scatin
        read(scatin_nodes_unit,*,iostat=istat) scatin_nodes(ii)
    end do
    
    close(scatin_nodes_unit)
    
    read(scatin_elements_unit,*,iostat=istat) m_scatin
    
    allocate(scatin_elements(m_scatin))
    
    do ii=1,m_scatin
        read(scatin_elements_unit,*,iostat=istat) scatin_elements(ii)
    end do
    
    close(scatin_elements_unit)
    
    read(scatb_elements_unit,*,iostat=istat) m_scatb
    
    allocate(scatb_elements(m_scatb))
    
    do ii=1,m_scatb
        read(scatb_elements_unit,*,iostat=istat) scatb_elements(ii)
    end do
    
    close(scatb_elements_unit)

    ! Huygens region
    
    read(huygb_nodes_unit,*,iostat=istat) n_huygb
    
    allocate(huygb_nodes(n_huygb))
    
    do ii=1,n_huygb
        read(huygb_nodes_unit,*,iostat=istat) huygb_nodes(ii)
    end do
    
    close(huygb_nodes_unit) 
    
    
    read(huygb_elements_unit,*,iostat=istat) m_huygb
    
    allocate(huygb_elements(m_huygb))
    
    do ii=1,m_huygb
        read(huygb_elements_unit,*,iostat=istat) huygb_elements(ii)
    end do
    
    close(huygb_elements_unit)  
    

                
end subroutine mesh_reader