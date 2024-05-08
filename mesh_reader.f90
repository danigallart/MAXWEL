subroutine mesh_reader
    use def_io
    use def_variables
    use def_vectors
    
    implicit none
    
    integer :: istat, ii
        
    double precision :: plasma_radius, fsdim, pmldim, huygdim, &
                        rpmlin, rpmlout, rhuyg
    double precision :: r
        
    ! Node's coordinates
    
    READ(nodes_unit,*,IOSTAT=istat) delh,NP
    
    allocate(coorx(NP), coory(NP))
    do ii=1,NP
        read(nodes_unit,*,iostat=istat) coorx(ii),coory(ii)
    end do
    
    CLOSE(nodes_unit)
        
    ! Conectivity matrix
    
    nodpel = 3
    
    read(elements_unit,*,iostat=istat) NE
    
    allocate(conn(NE,nodpel))
    DO ii=1,NE
        read(elements_unit,*,iostat=istat)  conn(ii,nodpel-2), conn(ii,nodpel-1), conn(ii,nodpel)
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
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate(coorx_mid(NE),coory_mid(NE))
    
    do ii=1,NE
    coorx_mid(ii) = sum(coorx(conn(ii,:)))/size(coorx(conn(ii,:)))
    coory_mid(ii) = sum(coory(conn(ii,:)))/size(coory(conn(ii,:)))
    enddo
    
    allocate(material(NE), boundary(NP))
    
    plasma_radius = r_scat * lambda0
    fsdim = lambda0/2
    pmldim = lambda0/2
    huygdim = delh
    
    rpmlin = plasma_radius + fsdim
    rpmlout = rpmlin + pmldim
    rhuyg = plasma_radius + huygdim
        
do ii=1,NP
    r = sqrt(coorx(ii)**2 + coory(ii)**2)
    if ((r < plasma_radius+1e-6) .and. (r > plasma_radius-1e-6)) then                       ! Boundary of scatterer
        boundary(ii) = 1
    elseif ((r < rhuyg+1e-6) .and. (r > rhuyg-1e-6)) then                                   ! Huygens surface
        boundary(ii) = 2
    elseif ((r < rpmlin+1e-6) .and. (r > rpmlin-1e-6)) then                                 ! Inner boundary of PML region
        boundary(ii) = 3
    elseif ((r < rpmlout+1e-6) .and. (r > rpmlout-1e-6)) then                               ! Outer boundary of PML region
        boundary(ii) = 4
    else                                                                                    ! No boundary node
        boundary(ii) = 0
    endif
enddo
    
do ii=1,NE
    r = sqrt(coorx_mid(ii)**2 + coory_mid(ii)**2)
    if (r < plasma_radius) then                                                             ! Scatterer element
        material(ii) = 1
    elseif ((r < rpmlout) .and. (r > rpmlin)) then                                          ! PML element
        material(ii) = 3
    else                                                                                    ! Vacuum element
        material(ii) = 2
    endif
enddo

        


end subroutine mesh_reader