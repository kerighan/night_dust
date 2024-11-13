from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
    name="night_dust",
    version="0.1.0",
    ext_modules=cythonize(
        Extension(
            name="night_dust.particles",
            sources=["night_dust/particles.pyx"],
            language_level="3",
        )
    ),
    packages=["night_dust"],
    zip_safe=False,
    install_requires=["cython"],
)
