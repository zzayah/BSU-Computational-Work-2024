from setuptools import setup, find_packages

setup(
    name="PVAssess",
    version="0.0.5",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'input_file_generation=pvskite.generation:create_input',
            'out_file_processing=pvskite.processing:process_output'
        ],
    
    },
    install_requires=[
        "ase",
        "pandas",
    ],
    author=["Zayah Cortright", "Stephen Lim"],
    author_email=["zayahcortright@gmail.com", "stephenlimzy2008@gmail.com"],
    description="",
    long_description=open("README.md").read(),
    url="",
)