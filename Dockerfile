
FROM python:3

COPY requirements.txt /


ADD ATAC-seqToolv1.5.py globalvariables.py /



RUN    pip install --no-cache-dir --upgrade pip && \
       pip install --no-cache-dir -r requirements.txt


CMD [ "python","./ATAC-seqToolv1.5.py" ]
