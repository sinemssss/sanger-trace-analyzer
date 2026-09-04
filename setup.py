from setuptools import setup, find_packages

setup(
    name="sanger-trace-analyzer",
    version="0.1.0",
    description="Sanger (.ab1) kromatogram pik analiz aracı - normal/heterozigot pozisyon tespiti",
    packages=find_packages(exclude=["tests", "examples"]),
    install_requires=[
        "pandas>=2.0",
        "matplotlib>=3.7",
    ],
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "sanger-analyze=sanger_trace_analyzer.main:main",
        ],
    },
)
