BUILD    := $(shell git describe --tags --always --dirty)
DIR := $(shell pwd)
TEST_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))test

# define image names
APP      := fhprs
REGISTRY := seglh

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(BUILD)
IMG_LATEST    := $(IMG):latest

.PHONY: push build tag test version cleanbuild

push: build
	docker push $(IMG_VERSIONED)
	docker push $(IMG_LATEST)

build: version
	docker buildx build --platform linux/amd64 -t $(IMG_VERSIONED) . || docker build -t $(IMG_VERSIONED) .
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)
	docker save $(IMG_VERSIONED) | gzip > $(DIR)/$(REGISTRY)-$(APP):$(BUILD).tar.gz

test: build
	echo "Incomplete VCF (positions missing):"
	docker run -it --rm -v $(TEST_DIR):/resources $(IMG_VERSIONED) /resources/incomplete.vcf
	echo "Complete VCF (all variant positions covered):"
	docker run -it --rm -v $(TEST_DIR):/resources $(IMG_VERSIONED) /resources/complete.vcf

cleanbuild:
	docker buildx build --platform linux/amd64 --no-cache -t $(IMG_VERSIONED) .