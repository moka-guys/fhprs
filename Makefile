BUILD    := $(shell git log -1 --pretty=%h)
TEST_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))test

# define image names
APP      := scidb
REGISTRY := fhprs

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(BUILD)
IMG_LATEST    := $(IMG):latest

.PHONY: push build tag test

.PHONY: build

push: build tag
	docker push $(IMG_VERSIONED)
	docker push $(IMG_LATEST)

build:
	docker build -t $(IMG_VERSIONED) .
	#docker buildx build --platform linux/amd64 -t $(IMG_VERSIONED) .

tag:
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)

test: build
	echo "Incomplete VCF (positions missing):"
	docker run -it --rm -v $(TEST_DIR):/resources $(IMG_VERSIONED) /resources/incomplete.vcf
	echo "Complete VCF (all variant positions covered):"
	docker run -it --rm -v $(TEST_DIR):/resources $(IMG_VERSIONED) /resources/complete.vcf

test_sample: build
	echo "Testing sample 17 from run NGS605:"
	docker run -it --rm -v $(TEST_DIR):/resources $(IMG_VERSIONED) /resources/NGS605_17_324382_LB_F_VCP1R134StG_Pan4821_S9_markdup_Haplotyper.g.vcf