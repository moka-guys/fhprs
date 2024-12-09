FROM python:3.6
RUN mkdir /code /sandbox /resources
WORKDIR /code
ADD ./fh.py /code
ADD ./test/tests_prs.py /code
ADD ./test/complete.vcf /code
ADD ./test/incomplete.vcf /code
ADD ./test/empty.vcf /code
ADD ./test/headers_only.vcf /code
ADD ./test/invalid.vcf /code
ADD requirements.txt /code
RUN pip3 install -r requirements.txt --retries=10
RUN chmod +x /code/fh.py
RUN python -m unittest -v tests_fhprs.py
ENTRYPOINT ["python","/code/fh.py"]

