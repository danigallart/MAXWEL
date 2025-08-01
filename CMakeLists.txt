cmake_minimum_required(VERSION 3.16)

# Nombre del proyecto
project(MAXWEL LANGUAGES Fortran)

# Activa advertencias de compilador (puedes quitar si molestan)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -fimplicit-none -ffree-line-length-none -Wno-error=line-truncation)
endif()

# Busca solo en la carpeta raíz del código
file(GLOB SRC_FILES CONFIGURE_DEPENDS
    *.f90
)

# Elimina duplicados por si acaso
list(REMOVE_DUPLICATES SRC_FILES)


# Crea el ejecutable
add_executable(maxwel ${SRC_FILES})


# ================================
# Configuración para MUMPS
# ================================

# Ruta base de instalación de MUMPS (ajusta según tu instalación)
set(MUMPS_DIR "/home/hernandr/MUMPS_5.6.2")

# Incluye cabeceras (mod files)
target_include_directories(maxwel PRIVATE
    ${MUMPS_DIR}/include
)

# Enlaza librerías de MUMPS y dependencias típicas
target_link_libraries(maxwel PRIVATE
    ${MUMPS_DIR}/lib/libdmumps.a
    ${MUMPS_DIR}/lib/libmumps_common.a
    ${MUMPS_DIR}/lib/libpord.a
    ${MUMPS_DIR}/lib/libmpiseq.a    # Quitar si usas versión MPI pura
    scalapack
    blas
    lapack
    mpi
)