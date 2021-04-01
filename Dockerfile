FROM python:3.7
RUN mkdir /code /sandbox /resources
WORKDIR /code
ADD ./fhprs /code
RUN pip3 install PyVCF==0.6.8
RUN chmod +x /code/fh.py
ENTRYPOINT ["python","/code/fh.py"]

