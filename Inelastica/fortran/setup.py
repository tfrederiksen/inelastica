def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    #from os.path import join as osp_join

    config = Configuration('fortran', parent_package, top_path)

    # Build the generic fortran helper routines
    # These do not depend on any external libraries
    all_info = get_info('ALL')
    sources = ['expansion_SE.f90',
               'readTSHS.f90',
               'removeUnitCellXij.f90',
               'setkpointhelper.f90']
    config.add_extension('F90helpers',
                         #sources=[osp_join('fortran', s) for s in sources],
                         sources=sources,
                         extra_info=all_info)

    # Build the lapack dependent routines.
    lapack_opt = get_info('lapack_opt')
    if not lapack_opt:
        raise NotFoundError('No LAPACK/BLAS resources found')
    sources = ['surfaceGreen.f90']
    config.add_extension('F90_lapack',
                         #sources=[osp_join('fortran', s) for s in sources],
                         sources=sources,
                         **lapack_opt)

    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
