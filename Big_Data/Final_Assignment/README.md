# Pneumonia-ID-via-X-Ray-images
Adrián Barreno Sánchez (adrian.barreno@alumnos.upm.es), Alberto González Delgado (alberto.gondelgado@alumnos.upm.es), Julian Elijah Politsch (julian.politsch@alumnos.upm.es), Angelo D'Angelo (angelo.dangelo@alumnos.upm.es)

This repository contains the resources of Big Data Engineering subject final project (MSc Computational Biology - Universidad Politécnica de Madrid).

All four members contributed equally, Julian did Neural Nework model, Alberto did Logistic Regression, Adrian did Random Forest Tree and Angelo did Support Vector Classification, and for that reason the evaluation should be equal for all members of the group.

# Aim

In this study, we aimed to utilize artificial intelligence (AI) to create a model for classifying chest X-ray images of individuals into two categories: pneumonia and healthy. With the increasing amount of medical imaging data being generated, manual interpretation by radiologists can become time-consuming and prone to errors. Our model utilizes a deep-learning algorithm, specifically convolutional neural networks (CNNs), to analyze the chest x-ray images and make predictions on the presence or absence of pneumonia. The performance of the model was evaluated using a dataset of chest X-ray images, and the results demonstrate the potential of AI in accurately classifying medical images and aiding in the diagnosis of pneumonia.

# Requirements

* [Scala Programing Languaje](https://www.scala-lang.org/)
* [sparkdl](https://gist.github.com/innat/b0ab252c954eb2a28a984774e3ee1f2d)
* [bigdl-dllib](https://sourceforge.net/projects/analytics-zoo/files/dllib-py-spark3/bigdl_dllib_spark3-0.14.0b20211107-py3-none-manylinux1_x86_64.whl)
* [Pandas](https://pandas.pydata.org/docs/getting_started/install.html)
* [Pillow](https://pypi.org/project/Pillow/)
* [NumPy](https://numpy.org/install/)
* [io](https://docs.python.org/es/3.9/library/io.html)
* [Tensorflow](https://pypi.org/project/tensorflow) 
* [PySpark](https://spark.apache.org/docs/latest/api/python/getting_started/install.html)

# Data

The data used in this project has been downloaded from the following [Kaggle link]( https://www.kaggle.com/paultimothymooney/chest-xray-pneumonia). This data belongs to a recent scientific publication from *Cell* journal, in which [Kermany et al. (2018)](https://www.cell.com/cell/fulltext/S0092-8674(18)30154-5) develoved a novel Artificial Intelligence algorithm to process medical images and provide an accurate diagnosis for different pathologies.

The dataset is organized into 3 folders (train, test, val) and contains subfolders for each image category (Pneumonia/Normal). There are 5,863 X-Ray images (JPEG) and 2 categories (Pneumonia/Normal).

Chest X-ray images (anterior-posterior) were selected from retrospective cohorts of pediatric patients of one to five years old from Guangzhou Women and Children’s Medical Center, Guangzhou. All chest X-ray imaging was performed as part of patients’ routine clinical care.
For the analysis of chest x-ray images, all chest radiographs were initially screened for quality control by removing all low quality or unreadable scans. The diagnoses for the images were then graded by two expert physicians before being cleared for training the AI system. In order to account for any grading errors, the evaluation set was also checked by a third expert.

# Analysis

The different steps of data analysis, processing, model building and evaluation are performed in the Jupyter Notebook [Pneumonia_Image_Classification.ipynb](./Pneumonia_Image_Classification.ipynb). The project pipeline can be summarized in 5 main steps:

1. We first load the data from Kaggle and pre-process it to obtain clean and tidy training and test datasets.
2. Images are transformed into numeric arrays that represent the image features and are easy interpretable by the machine learning models.
3. Principal Component Analysis to reduce dymensionality and time consumption of the problem.
4. Train the different models and perform hyperparameter tuning via cross-validation.
5. Select the best performing model and use it to predict the status of the test dataset images.