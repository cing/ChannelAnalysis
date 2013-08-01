from distutils.core import setup

setup(
    name='ChannelAnalysis',
    version='0.1.0',
    author='Christopher Ing',
    author_email='ing.chris@gmail.com',
    packages=['ChannelAnalysis','ChannelAnalysis.RotamerAnalysis',
              'ChannelAnalysis.PoreAnalysis','ChannelAnalysis.CoordAnalysis'],
    url='https://github.com/cing/ChannelAnalysis/',
    license='LICENSE',
    description='An analysis pipeline for studying ion channels.',
    long_description=open('README.markdown').read(),
    install_requires=[
        "numpy",
        "scipy",
    ],
)
