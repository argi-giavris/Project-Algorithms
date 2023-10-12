import random
import sys
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False

from keras.models import load_model
import tensorflow
import pandas as pd
import numpy
from tensorflow import keras
from keras import layers, optimizers, losses, metrics
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout
import math
from sklearn.preprocessing import MinMaxScaler
import sys
import matplotlib.pyplot as plt

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

    number_of_epochs = 2
    layer_size = 20
    lookback = 1
    number_of_layers = 4
    batch_size = 1

    df = pd.read_csv(data_path,sep = '\t', header = None)
    print('Number of rows and columns:', df.shape)
    print('Number of rows and columns:', df.shape)

    y = 0
    while(y < number_of_timeseries):
        stock = random.randrange(df.shape[0])
        print(stock)
        data = df.iloc[[stock]]
        print(df[0][y])


        y = y+1
        data = data.T
        training_rows = math.floor(data.shape[0] * 0.8)

        training_set = data.iloc[1:training_rows, :].values
        test_set = data.iloc[training_rows:, :].values


        # Feature Scaling
        sc = MinMaxScaler(feature_range=(0, 1))
        training_set_scaled = sc.fit_transform(training_set)
        # Creating a data structure with 60 time-steps and 1 output
        X_train = []
        y_train = []
        for i in range(lookback, training_rows-1):
            X_train.append(training_set_scaled[i - lookback:i, 0])
            y_train.append(training_set_scaled[i, 0])
        X_train, y_train = numpy.array(X_train), numpy.array(y_train)
        X_train = numpy.reshape(X_train, (X_train.shape[0], X_train.shape[1], 1))
        # (740, 60, 1)

        model = Sequential()
        # Adding the first LSTM layer and some Dropout regularisation
        model.add(LSTM(units=layer_size, return_sequences=True, input_shape=(X_train.shape[1], 1)))
        model.add(Dropout(0.2))
        # Adding a second LSTM layer and some Dropout regularisation
        x = 0
        while number_of_layers - 3 > x:
            model.add(LSTM(units=layer_size, return_sequences=True))
            model.add(Dropout(0.2))
            x = x+1

        # Adding a fourth LSTM layer and some Dropout regularisation
        model.add(LSTM(units=layer_size))
        model.add(Dropout(0.2))

        #model.add(LSTM(4, input_shape=(1, look_back)))

        model.add(Dense(1))
        model.compile(loss='mean_squared_error', optimizer='adam')
        model.fit(X_train, y_train, epochs=number_of_epochs, batch_size=batch_size, verbose=2)

        # Getting the predicted stock price of 2017
        dataset_train = data.iloc[:training_rows, :]
        dataset_test = data.iloc[training_rows:, :]
        dataset_total = pd.concat((dataset_train, dataset_test), axis=0)
        inputs = dataset_total[len(dataset_total) - len(dataset_test) - 60:].values
        inputs = inputs.reshape(-1, 1)
        inputs = sc.transform(inputs)
        X_test = []
        test_rows = test_set.shape[0]
        for i in range(lookback, test_set.shape[0]+lookback):
            X_test.append(inputs[i - lookback:i, 0])
        X_test = numpy.array(X_test)
        X_test = numpy.reshape(X_test, (X_test.shape[0], X_test.shape[1], 1))
        print(X_test.shape)
        # (459, 60, 1)

        predicted_stock_price = model.predict(X_test)
        predicted_stock_price = sc.inverse_transform(predicted_stock_price)

        # Visualising the results
        # shift test predictions for plotting
        # Visualising the results
        plt.plot( dataset_test.values, color = 'red', label = 'Real')
        plt.plot( predicted_stock_price, color = 'blue', label = 'Predicted')
        plt.xticks(numpy.arange(0, test_rows, 50))
        plt.title(df[0][stock])
        plt.xlabel('Time')
        plt.ylabel('TESLA Stock Price')
        plt.legend()
        plt.show()

    y = 0
    sc = MinMaxScaler(feature_range=(0, 1))
    while(y < df.shape[0]):
        stock = y
        data = df.iloc[[stock]]


        y = y+1
        data = data.T
        training_rows = math.floor(data.shape[0] * 0.8)

        training_set = data.iloc[1:training_rows, :].values
        test_set = data.iloc[training_rows:, :].values



        # Feature Scaling
        training_set_scaled = sc.fit_transform(training_set)
        # Creating a data structure with 60 time-steps and 1 output
        if y != 1:
            previousX_train = X_train
            previousy_train = y_train
        X_train = []
        y_train = []
        for i in range(lookback, training_rows-1):
            X_train.append(training_set_scaled[i - lookback:i, 0])
            y_train.append(training_set_scaled[i, 0])
        X_train, y_train = numpy.array(X_train), numpy.array(y_train)
        X_train = numpy.reshape(X_train, (X_train.shape[0], X_train.shape[1], 1))
        if y != 1:
            X_train = numpy.concatenate((X_train, previousX_train), axis=0)
            y_train = numpy.concatenate((y_train, previousy_train), axis=0)
        # (740, 60, 1)

    model = Sequential()
    # Adding the first LSTM layer and some Dropout regularisation
    model.add(LSTM(units=layer_size, return_sequences=True, input_shape=(X_train.shape[1], 1)))
    model.add(Dropout(0.2))
    # Adding a second LSTM layer and some Dropout regularisation
    x = 0
    while number_of_layers - 3 > x:
        model.add(LSTM(units=layer_size, return_sequences=True))
        model.add(Dropout(0.2))
        x = x+1

    # Adding a fourth LSTM layer and some Dropout regularisation
    model.add(LSTM(units=layer_size))
    model.add(Dropout(0.2))

    #model.add(LSTM(4, input_shape=(1, look_back)))

    model.add(Dense(1))
    model.compile(loss='mean_squared_error', optimizer='adam')
    model.fit(X_train, y_train, epochs=number_of_epochs, batch_size=batch_size, verbose=2)

    y = 0
    while (y < number_of_timeseries):
        y = y+1
        stock = random.randrange(df.shape[0])
        print(stock)
        data = df.iloc[[stock]]
        data = data.T
        print(df[0][stock])

        # Getting the predicted stock price of 2017
        dataset_train = data.iloc[:training_rows, :]
        dataset_test = data.iloc[training_rows:, :]
        dataset_total = pd.concat((dataset_train, dataset_test), axis=0)
        inputs = dataset_total[len(dataset_total) - len(dataset_test) - 60:].values
        inputs = inputs.reshape(-1, 1)
        inputs = sc.transform(inputs)
        X_test = []
        test_rows = test_set.shape[0]
        for i in range(lookback, test_set.shape[0]+lookback):
            X_test.append(inputs[i - lookback:i, 0])
        X_test = numpy.array(X_test)
        X_test = numpy.reshape(X_test, (X_test.shape[0], X_test.shape[1], 1))
        print(X_test.shape)
        # (459, 60, 1)

        predicted_stock_price = model.predict(X_test)
        predicted_stock_price = sc.inverse_transform(predicted_stock_price)

        # Visualising the results
        # shift test predictions for plotting
        # Visualising the results
        plt.plot( dataset_test.values, color = 'red', label = 'Real')
        plt.plot( predicted_stock_price, color = 'blue', label = 'Predicted')
        plt.xticks(numpy.arange(0, test_rows, 50))
        plt.title(df[0][stock])
        plt.xlabel('Time')
        plt.ylabel('TESLA Stock Price')
        plt.legend()
        plt.show()
