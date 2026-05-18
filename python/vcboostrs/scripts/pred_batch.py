import tensorflow as tf
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Layer, Flatten, Dense, Bidirectional, Dropout, LSTM, \
    LayerNormalization, MultiHeadAttention, Concatenate
from tensorflow.keras.models import Model
from tensorflow.keras import backend as K
import sys
import time
import logging
import os
from tensorflow.keras.layers import Input, Conv2D, BatchNormalization, Activation, Add
import pandas as pd
import numpy as np

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

log_path = "predict.log"
if os.path.exists(log_path):
    index = 1
    while True:
        new_path = f"predict_{index}.log"
        if not os.path.exists(new_path):
            log_path = new_path
            break
        index += 1

logger = logging.getLogger('predict_logger')
logger.setLevel(logging.INFO)
file_handler = logging.FileHandler(log_path)
formatter = logging.Formatter('%(asctime)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


def residual_block(input_tensor, filters):
    x = Conv2D(filters, kernel_size=(3, 3), padding='same')(input_tensor)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x = Conv2D(filters, kernel_size=(3, 3), padding='same')(x)
    x = BatchNormalization()(x)
    x = Add()([x, input_tensor])
    x = Activation('relu')(x)
    return x


class DifferenceLayer(Layer):
    def __init__(self, **kwargs):
        super(DifferenceLayer, self).__init__(**kwargs)

    def call(self, inputs):
        if inputs.shape[-1] != 4:
            raise ValueError("DifferenceLayer expects the last dimension to be 4.")
        diff1 = K.abs(inputs[:, :, :, 0] - inputs[:, :, :, 2])
        diff2 = K.abs(inputs[:, :, :, 1] - inputs[:, :, :, 3])
        diff3 = K.abs(inputs[:, :, :, 0] - inputs[:, :, :, 1])
        diff4 = K.abs(inputs[:, :, :, 2] - inputs[:, :, :, 3])
        merged_diff = K.stack([diff1, diff2, diff3, diff4], axis=-1)
        return merged_diff


def residual_convnet_model(input_shape, num_classes):
    inputs = Input(shape=input_shape)
    x = Conv2D(128, kernel_size=(7, 5), padding='same')(inputs)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x = make_basic_block_layer(filter_num=128, blocks=1, stride=1)(x)
    x = Conv2D(128, kernel_size=(7, 5), strides=(1, 1), padding="same")(x)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x = make_basic_block_layer(filter_num=128, blocks=1, stride=1)(x)
    x = Flatten()(x)
    x = tf.keras.layers.Reshape((72, -1))(x)
    x = Bidirectional(LSTM(128, return_sequences=True))(x)
    x = MultiHeadAttention(key_dim=128, num_heads=8)(x, x)
    x = Bidirectional(LSTM(128, return_sequences=True))(x)
    x = Bidirectional(LSTM(128, return_sequences=True))(x)
    x = LSTM(128)(x)
    x = Dense(units=256, activation='relu')(x)
    x = Dropout(0.5)(x)
    outputs = Dense(units=num_classes, activation='softmax')(x)
    model = Model(inputs=inputs, outputs=outputs)
    return model


L2_regularizers = tf.keras.regularizers.l2(1e-7)


class BasicBlock(tf.keras.layers.Layer):
    def __init__(self, filter_num, stride=1):
        super(BasicBlock, self).__init__()
        self.filter_num = filter_num
        self.stride = stride
        self.conv1 = tf.keras.layers.Conv2D(filters=filter_num, kernel_size=(7, 5), strides=stride, padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn1 = tf.keras.layers.BatchNormalization()
        self.conv2 = tf.keras.layers.Conv2D(filters=filter_num, kernel_size=(7, 5), strides=1, padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn2 = tf.keras.layers.BatchNormalization()
        if stride != 1:
            self.downsample = tf.keras.Sequential()
            self.downsample.add(tf.keras.layers.Conv2D(filters=filter_num, kernel_size=(7, 5), strides=stride,
                                                       kernel_regularizer=L2_regularizers))
            self.downsample.add(tf.keras.layers.BatchNormalization())
        else:
            self.downsample = lambda x: x

    def call(self, inputs):
        residual = self.downsample(inputs)
        x = self.conv1(inputs)
        x = self.bn1(x, )
        x = tf.nn.relu(x)
        x = self.conv2(x)
        x = self.bn2(x, )
        output = tf.nn.relu(tf.keras.layers.add([residual, x]))
        return output

    def get_config(self):
        config = super().get_config()
        config.update({'filter_num': self.filter_num, 'stride': self.stride})
        return config


def make_basic_block_layer(filter_num, blocks, stride=1):
    res_block = tf.keras.Sequential()
    res_block.add(BasicBlock(filter_num, stride=stride))
    for _ in range(1, blocks):
        res_block.add(BasicBlock(filter_num, stride=1))
    return res_block


pred = sys.argv[1]
input_shape = (72, 5, 4)
num_classes = 2

logger.info(f"pred_batch.py started")
logger.info(f"Input file: {pred}")

t_start = time.time()

model = residual_convnet_model(input_shape, num_classes)
path = sys.argv[2]
logger.info(f"Model weights: {path}")

model.load_weights(path)
logger.info("Model weights loaded")

t_load = time.time()
logger.info(f"Model loading time: {t_load - t_start:.2f}s")

data = pd.read_csv(pred, sep=' ', header=None)
num_rows = len(data)
logger.info(f"Data loaded: {num_rows} rows")

chr_col = data.iloc[:, 0]
pos = data.iloc[:, 1]
data_values = data.iloc[:, 2:].values.astype(np.int32)

if data_values.shape[1] > 72 * 5 * 4:
    data_values = data_values[:, :72 * 5 * 4]

data_reshaped = data_values.reshape(-1, 72, 5, 4)

batch = 256
logger.info(f"Predicting {num_rows} samples with batch_size={batch}...")

t_pred_start = time.time()
pred_labels_one = model.predict(data_reshaped, batch_size=batch)
t_pred_end = time.time()
logger.info(f"Prediction time: {t_pred_end - t_pred_start:.2f}s")
logger.info(f"Prediction speed: {num_rows / (t_pred_end - t_pred_start):.1f} samples/s")

chr_out = chr_col.to_numpy().reshape(-1, 1)
pos_out = pos.to_numpy().reshape(-1, 1)
output_file_path = pred + '.txt'

out_array = np.concatenate([chr_out, pos_out, pred_labels_one], axis=1)
df = pd.DataFrame(out_array)
df.to_csv(output_file_path, sep=' ', index=False, header=False)

t_end = time.time()
logger.info(f"Total time: {t_end - t_start:.2f}s")
logger.info(f"Output saved to: {output_file_path}")
