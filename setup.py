from setuptools import setup, find_packages

setup(
    name='DrawTreesETE',
    version='0.1.13',
    author='Veronika Berman',
    author_email='nika.berman@gmail.com',
    packages = find_packages(),  # include all packages under src
    # package_dir = {'':'src'},   # tell distutils packages are under src
    #url='http://pypi.python.org/pypi/TowelStuff/',
    #license='LICENSE.txt',
    description='Useful towel-related stuff.',
    #long_description=open('README.txt').read(),
    install_requires=[
        "ete2 == 2.3.10.1",
    ],
    dependency_links=[
        "git+https://github.com/NikiB/ete.git@476722ec18338eb46037b06334e3b16c81bface2#egg=ete2-2.3.10.1",
    ],
    entry_points={
        'console_scripts': [
            'draw_tree=tree_viewer.draw_tree:main',
        ]
    }
)
