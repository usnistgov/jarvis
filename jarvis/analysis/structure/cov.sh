#!/bin/bash
pytest 
coverage run -m pytest
coverage report -m
