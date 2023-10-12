from keras.layers import Input, Dense, Conv1D, MaxPooling1D, UpSampling1D, BatchNormalization, LSTM, RepeatVector
from keras.models import Model
from keras.models import model_from_json
from keras import regularizers
import datetime
import time
import requests as req
import json
import pandas as pd
import pickle
import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import math
import random


def mkdate(ts):
    return datetime.datetime.fromtimestamp(
        int(ts)
    ).strftime('%Y-%m-%d')


def plot_examples(stock_input, stock_decoded):
    n = 10
    plt.figure(figsize=(20, 4))
    for i, idx in enumerate(list(np.arange(0, test_samples, 200))):
        # display original
        ax = plt.subplot(2, n, i + 1)
        if i == 0:
            ax.set_ylabel("Input", fontweight=600)
        else:
            ax.get_yaxis().set_visible(False)
        plt.plot(stock_input[idx])
        ax.get_xaxis().set_visible(False)

        # display reconstruction
        ax = plt.subplot(2, n, i + 1 + n)
        if i == 0:
            ax.set_ylabel("Output", fontweight=600)
        else:
            ax.get_yaxis().set_visible(False)
        plt.plot(stock_decoded[idx])
        ax.get_xaxis().set_visible(False)


def plot_history(history):
    plt.figure(figsize=(15, 5))
    ax = plt.subplot(1, 2, 1)
    plt.plot(history.history["loss"])
    plt.title("Train loss")
    ax = plt.subplot(1, 2, 2)
    plt.plot(history.history["val_loss"])
    plt.title("Test loss")

if __name__ == "__main__":
    arguments = len(sys.argv)
    print(arguments)
    j = 0
    for i in sys.argv:
        j = j + 1
        if(i == '-d'):
            data_path = sys.argv[j]
        if (i == '-n'):
            number_of_timeseries = int(sys.argv[j])
        if (i == '-mae'):
            mean_absolute_error = int(sys.argv[j])

    number_of_epochs = 5
    layer_size = 128
    number_of_layers = 4
    batch_size = 128
    TIME_STEPS = 30
    window_length = 10

    df = pd.read_csv(data_path,sep = '\t', header = None)
    print('Number of rows and columns:', df.shape)
    print('Number of rows and columns:', df.shape)
    data = df.T

    training_rows = math.floor(data.shape[0] * 0.8)

    y = 0
    while(y < df.shape[0]):
        stock = y
        data = df.iloc[[stock]]
        data = data.T


        y = y+1
        training_rows = math.floor(data.shape[0] * 0.8)

        #training_set = data.iloc[1:training_rows, :].values
        #test_set = data.iloc[training_rows:, :].values

        data = pd.DataFrame(np.array(data)[1:, 0], columns=['price'])
        data['pct_change'] = data.price.pct_change()

        data['log_ret'] = np.log(data['price'].astype('float')) - np.log(data['price'].astype('float').shift(1))

        scaler = MinMaxScaler()
        if y != 1:
          previousX_train = X_train
          previousy_train = y_train
        X_train = []
        y_train = []
        x_train_nonscaled = np.array([data['log_ret'].values[i - window_length:i].reshape(-1, 1) for i in
                                      tqdm(range(window_length + 1, len(data['log_ret'])))])
        x_train = np.array([scaler.fit_transform(data['log_ret'].values[i - window_length:i].reshape(-1, 1)) for i in
                            tqdm(range(window_length + 1, len(data['log_ret'])))])

        x_test = x_train[training_rows:]
        x_train = x_train[1:training_rows]

        x_train = x_train.astype('float32')
        x_test = x_test.astype('float32')

        if y != 1:
            X_train = numpy.concatenate((X_train, previousX_train), axis=0)
            y_train = numpy.concatenate((y_train, previousy_train), axis=0)


    input_window = Input(shape=(window_length,1))
    x = Conv1D(16, 3, activation="relu", padding="same")(input_window) # 10 dims
    #x = BatchNormalization()(x)
    x = MaxPooling1D(2, padding="same")(x) # 5 dims
    x = Conv1D(1, 3, activation="relu", padding="same")(x) # 5 dims
    #x = BatchNormalization()(x)
    encoded = MaxPooling1D(2, padding="same")(x) # 3 dims

    encoder = Model(input_window, encoded)

    # 3 dimensions in the encoded layer

    x = Conv1D(1, 3, activation="relu", padding="same")(encoded) # 3 dims
    #x = BatchNormalization()(x)
    x = UpSampling1D(2)(x) # 6 dims
    x = Conv1D(16, 2, activation='relu')(x) # 5 dims
    #x = BatchNormalization()(x)
    x = UpSampling1D(2)(x) # 10 dims
    decoded = Conv1D(1, 3, activation='sigmoid', padding='same')(x) # 10 dims
    autoencoder = Model(input_window, decoded)
    autoencoder.summary()

    autoencoder.compile(optimizer='adam', loss='binary_crossentropy')
    history = autoencoder.fit(x_train, x_train,
                    epochs=number_of_epochs,
                    batch_size=batch_size,
                    shuffle=True,
                    validation_data=(x_test, x_test))

    decoded_stocks = autoencoder.predict(x_test)
