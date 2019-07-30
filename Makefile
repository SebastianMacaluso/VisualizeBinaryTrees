# Makefile for VisualizeBinaryTrees
SHELL := /bin/bash

# You can set these variables from the commandline.
VERSION=$(shell python setup.py --version)

./dist/VisualizeBinaryTrees-${VERSION}-py3-none-any.whl:
	python ./setup.py sdist bdist_wheel

install: ./dist/VisualizeBinaryTrees-${VERSION}-py3-none-any.whl # pip install
	pip install --upgrade ./dist/VisualizeBinaryTrees-${VERSION}-py3-none-any.whl


%: Makefile