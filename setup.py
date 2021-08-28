from setuptools import setup, find_namespace_packages

setup(
    name="mirna_profiling",
    version="1.0.0",
    description="miRNA expression profiling",
    packages=find_namespace_packages(),
    install_requires=[
        "pytest",
        "pandas==1.2.5",
        "matplotlib",
        "pysam",
        "pyyaml",
        "PyMySQL",
        # "mysql-connector-python",
        "scipy",
        "statsmodels",
        "zc.buildout",
    ],
    package_data={'mirna_profiling': ['configuration/annotation_config.yaml', 'annotation/resource/*']},
    dependency_links=["https://pypi.bcgsc.ca/gsc/packages/"],
    url="",
    license="",
    author="ahe",
    python_requires=">=3.8",
    zip_safe=False,
    author_email="ahe@bcgsc.ca",
    entry_points={
        "console_scripts": [
            "annotate_bam=mirna_profiling.annotation.annotate_bam:main",
            "alignment_stat=mirna_profiling.library_stat.alignment_stat:main",
            "expression_matrix=mirna_profiling.library_stat.expression_matrix:main",
            "project_summary=mirna_profiling.library_stat.project_summary:main",
            "tcga=mirna_profiling.library_stat.tcga:main",
            "graph_stat=mirna_profiling.library_stat.graph_stat:main",
        ]
    },
)
