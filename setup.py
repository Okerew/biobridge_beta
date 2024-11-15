from setuptools import setup, Extension, find_packages
import pybind11
from glob import glob
import os


def get_all_src_files():
    return [file for file in glob('biobridge/**/*.py', recursive=True)]


ext_modules = [
    Extension(
        'biobridge.visualizer',
        ['biobridge/cpp/visualizer.cpp'],
        include_dirs=[
            pybind11.get_include(),
            '/opt/homebrew/Cellar/sfml/2.6.1/include'
        ],
        library_dirs=[
            '/opt/homebrew/Cellar/sfml/2.6.1/lib'
        ],
        libraries=['sfml-graphics', 'sfml-window', 'sfml-system'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        'biobridge.neural_network',
        ['biobridge/cpp/neural_network.cpp'],
        include_dirs=[pybind11.get_include(),
                      '/opt/homebrew/Cellar/nlohmann-json/3.11.3/include'
                      ],
        library_dirs=['/opt/homebrew/Cellar/nlohmann-json/3.11.3/lib'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        name='biobridge.operations_simulator',
        sources=['biobridge/cpp/operations_simulator.cpp'],
        include_dirs=[pybind11.get_include(),
                    '/opt/homebrew/Cellar/sfml/2.6.1/include'
                ],
        library_dirs=[
            '/opt/homebrew/Cellar/sfml/2.6.1/lib'
        ],
        libraries=['sfml-graphics', 'sfml-window', 'sfml-system'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        name="biobridge.metabolic_network",
        sources=["biobridge/cpp/metabolic_network.cpp"],
        include_dirs=[pybind11.get_include(),
                     '/opt/homebrew/Cellar/nlohmann-json/3.11.3/include',
                    '/opt/homebrew/Cellar/boost/1.86.0/include'
                    ],
        library_dirs=['/opt/homebrew/Cellar/nlohmann-json/3.11.3/lib',
                      '/opt/homebrew/Cellar/boost/1.86.0/lib'],
        libraries=['boost_filesystem', 'boost_system'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        name="biobridge.signaling_network",
        sources=["biobridge/cpp/signaling_network.cpp"],
        include_dirs=[pybind11.get_include(),
                     '/opt/homebrew/Cellar/nlohmann-json/3.11.3/include',
                    ],
        library_dirs=['/opt/homebrew/Cellar/nlohmann-json/3.11.3/lib'],
        libraries=['boost_filesystem', 'boost_system'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        name="biobridge.individuals",
        sources=["biobridge/cpp/consciousness_simulator.cpp", "biobridge/cpp/society_simulator.cpp"],
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        name="biobridge.immune_system",
        sources=["biobridge/cpp/immune_system.cpp"],
        include_dirs=[pybind11.get_include(), '/opt/homebrew/Cellar/sfml/2.6.1/include'],
        library_dirs=['/opt/homebrew/Cellar/sfml/2.6.1/lib'],
        libraries=['sfml-graphics', 'sfml-window', 'sfml-system'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        name="biobridge.embryo_simulator",
        sources=["biobridge/cpp/embryo_simulator.cpp"],
        include_dirs=[pybind11.get_include(), '/opt/homebrew/Cellar/sfml/2.6.1/include'],
        library_dirs=['/opt/homebrew/Cellar/sfml/2.6.1/lib'],
        libraries=['sfml-graphics', 'sfml-window', 'sfml-system'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
    Extension(
        name="biobridge.infection_simulator",
        sources=["biobridge/cpp/infection_simulator.cpp"],
        include_dirs=[pybind11.get_include(), '/opt/homebrew/Cellar/sfml/2.6.1/include'],
        library_dirs=['/opt/homebrew/Cellar/sfml/2.6.1/lib'],
        libraries=['sfml-graphics', 'sfml-window', 'sfml-system'],
        language='c++',
        extra_compile_args=['-std=c++11'],
    ),
]

os.system("python create_config.py")

setup(
    name='biobridge',
    version='0.0.1',
    author='Okerew',
    author_email='okerewgroup@proton.me',
    description='A visualizer for the cell environment',
    packages=find_packages(),
    package_data={'': ['*.txt', '*.json']},
    include_package_data=True,
    ext_modules=ext_modules,
    install_requires=['pybind11'],
    py_modules=[os.path.splitext(os.path.basename(path))[0] for path in get_all_src_files()],
)
