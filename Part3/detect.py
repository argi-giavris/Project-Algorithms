import random
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False

import pandas as pd
import numpy
from tensorflow import keras
import seaborn as sns
import math
import sys
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

def create_dataset(X, y, time_steps=1):
    Xs, ys = [], []
    for i in range(len(X) - time_steps):
        v = X.iloc[i:(i + time_steps)].values
        Xs.append(v)
        ys.append(y.iloc[i + time_steps])
    return numpy.array(Xs), numpy.array(ys)

if __name__ == "__main__":
    arguments = len(sys.argv)
    print(arguments)
    mean_absolute_error = 0.65
    j = 0
    for i in sys.argv:
        j = j + 1
        if(i == '-d'):
            data_path = sys.argv[j]
        if (i == '-n'):
            number_of_timeseries = int(sys.argv[j])
        if (i == '-mae'):
            mean_absolute_error = int(sys.argv[j])

    number_of_epochs = 10
    layer_size = 64
    number_of_layers = 4
    batch_size = 128
    TIME_STEPS = 30

    df = pd.read_csv(data_path,sep = '\t', header = None)
    print('Number of rows and columns:', df.shape)
    print('Number of rows and columns:', df.shape)
    data = df.T

    training_rows = math.floor(data.shape[0] * 0.8)

    y = 0
    date = [*range(0,data.shape[0])]
    print(date)
    scaler = StandardScaler()
    while(y < df.shape[0]):
        stock = y
        data = df.iloc[[stock]]
        data = data.T
        data['date'] = date
        y = y + 1
        training_set = data.iloc[1:training_rows]
        training_set.rename(columns={ y-1 : 'Close'},inplace = True)
        '''plt.plot(training_set['Close'], label='Close')
        plt.legend();
        #plt.show()'''

        #test_set = data.iloc[training_rows:]
        #scaler = StandardScaler()
        scaler = scaler.fit(training_set[['Close']])
        training_set['Close'] = scaler.transform(training_set[['Close']])
        #test_set = scaler.transform(test_set)

        if y != 1:
            previousX_train = X_train
            previousy_train = y_train
        X_train = []
        y_train = []
        # reshape to [samples, time_steps, n_features]
        X_train, y_train = create_dataset(training_set[['Close']],training_set.Close,TIME_STEPS)
        #X_test, y_test = create_dataset(test_set,range(training_rows,training_rows+test_set.shape[0]),TIME_STEPS)
        if y != 1:
            X_train = numpy.concatenate((X_train, previousX_train), axis=0)
            y_train = numpy.concatenate((y_train, previousy_train), axis=0)


    print(X_train.shape)
    model = keras.Sequential()
    model.add(keras.layers.LSTM(units=layer_size,input_shape=(X_train.shape[1], X_train.shape[2])))
    model.add(keras.layers.Dropout(rate=0.2))
    model.add(keras.layers.RepeatVector(n=X_train.shape[1]))
    model.add(keras.layers.LSTM(units=layer_size, return_sequences=True))
    model.add(keras.layers.Dropout(rate=0.2))
    model.add(keras.layers.TimeDistributed(keras.layers.Dense(units=X_train.shape[2])))
    model.compile(loss='mae', optimizer='adam')

    history = model.fit(X_train, y_train,epochs=number_of_epochs,batch_size=batch_size,validation_split=0.1,shuffle=False)

    '''plt.plot(history.history['loss'], label='train')
    plt.plot(history.history['val_loss'], label='test')
    plt.legend();
    plt.show()
    print("kati")'''
    y = 0
    while (y < number_of_timeseries):
        y = y+1
        stock = random.randrange(df.shape[0])
        stock = y-1
        data = df.iloc[[stock]]
        data = data.T
        data['date'] = date
        data.rename(columns={stock: 'Close'}, inplace=True)
        #training_set = data.iloc[1:training_rows]
        #print(training_set)
        test_set = data.iloc[training_rows:]
        #print(test_set)
        test_set['Close'] = scaler.transform(test_set[['Close']])
        X_test, y_test = create_dataset(test_set[['Close']], test_set.Close, TIME_STEPS)
        #scaler = scaler.fit(training_set[['Close']])
        #training_set['Close'] = scaler.transform(training_set[['Close']])
        #X_train, y_train = create_dataset(training_set[['Close']], training_set.Close, TIME_STEPS)


    #X_train_pred = model.predict(X_train)
    #train_mae_loss = numpy.mean(numpy.abs(X_train_pred - X_train), axis=1)

        THRESHOLD = mean_absolute_error

        X_test_pred = model.predict(X_test)
        test_mae_loss = numpy.mean(numpy.abs(X_test_pred - X_test), axis=1)
        test_score_df = pd.DataFrame(index=test_set[TIME_STEPS:].index)
        test_score_df['loss'] = test_mae_loss
        test_score_df['threshold'] = THRESHOLD
        test_score_df['anomaly'] = test_score_df.loss > test_score_df.threshold
        test_score_df['close'] = test_set[TIME_STEPS:].Close

        anomalies = test_score_df[test_score_df.anomaly == True]
        anomalies.head()

        plt.plot(
            test_set[TIME_STEPS:].index,
            scaler.inverse_transform(test_set[TIME_STEPS:].Close),
            label='close price'
        );

        sns.scatterplot(
            anomalies.index,
            scaler.inverse_transform(anomalies.close),
            color=sns.color_palette()[3],
            s=52,
            label='anomaly'
        )
        plt.xticks(rotation=25)
        plt.legend();
        plt.show()
