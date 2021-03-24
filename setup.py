from setuptools import setup

setup(name = "phononTools",
      version = "0.0.1",
      author = "Harrison LaBollita",
      author_email = "hlabolli@asu.edu",
      description = "a small python package for analyzing phonon calculations",
      packages = ["phononTools", "phononTools.base"],
      install_requires =["numpy", "matplotlib"],
      scripts = ["scripts/phononTools"])

