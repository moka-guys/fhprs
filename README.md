# FH-PRS

Simple standalone class to calculate FH polygenic risk score accoding to Talmund et al. 2013.
Code adopted from Bristol and modernised for python3 and pyVCF. Uses build 38 genomic coordinates. 

## Usage
`python3 fh.py file.vcf`

## Docker image build
`make build`

To build the docker image from scratch, using no cached layers:
`make cleanbuild` 

## Docker image upload
To upload the built docker image to DockerHub:
`make push`

## Test
This runs against the two files with all or partially missing positions
`make test`
