from setuptools import setup, find_packages

setup(
    name="pvskite-assess",
    version="0.0.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'input_file_generation=pvskite.generation:create_input'
            'out_file_processing=pvskite.processing:process_output'
        ],
    
    },
    install_requires=[
        "ase",
        "pandas",
    ],
    author="Zayah Cortright",
    author_email="zayahcortright@gmail.com",
    description="",
    long_description=open("README.md").read(),
    url="",
)