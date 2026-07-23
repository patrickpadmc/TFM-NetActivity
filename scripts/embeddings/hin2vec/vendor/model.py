"""
model.py -- vendorizado de e-yi/hin2vec_pytorch (MIT license), sin
modificaciones. No necesita parche: HIN2vec y NSTrainSet consumen
(start_node, end_node, path) como ids enteros opacos, sin importarles
que representa cada uno. La logica de que representa "path" vive en
el walker.py parchado (ver ese archivo).
"""
import torch
import numpy as np
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader


def binary_reg(x: torch.Tensor):
    # forward: f(x) = (x>=0)
    # backward: f(x) = sigmoid
    a = torch.sigmoid(x)
    b = a.detach()
    c = (x.detach() >= 0).float()
    return a - b + c


class HIN2vec(nn.Module):
    def __init__(self, node_size, path_size, embed_dim, sigmoid_reg=False, r=True):
        super().__init__()
        self.reg = torch.sigmoid if sigmoid_reg else binary_reg
        self.__initialize_model(node_size, path_size, embed_dim, r)

    def __initialize_model(self, node_size, path_size, embed_dim, r):
        self.start_embeds = nn.Embedding(node_size, embed_dim)
        self.end_embeds = self.start_embeds if r else nn.Embedding(node_size, embed_dim)
        self.path_embeds = nn.Embedding(path_size, embed_dim)

    def forward(self, start_node: torch.LongTensor, end_node: torch.LongTensor, path: torch.LongTensor):
        s = self.start_embeds(start_node)  # (batch_size, embed_size)
        e = self.end_embeds(end_node)
        p = self.path_embeds(path)
        p = self.reg(p)

        agg = torch.mul(s, e)
        agg = torch.mul(agg, p)

        output = torch.sigmoid(torch.sum(agg, axis=1))
        return output


def train(log_interval, model, device, train_loader: DataLoader, optimizer, loss_function, epoch):
    model.train()
    for idx, (data, target) in enumerate(train_loader):
        data, target = data.to(device), target.to(device)
        optimizer.zero_grad()
        output = model(data[:, 0], data[:, 1], data[:, 2])
        loss = loss_function(output.view(-1), target)
        loss.backward()
        optimizer.step()

        if idx % log_interval == 0:
            print(f'\rTrain Epoch: {epoch} '
                  f'[{idx * len(data)}/{len(train_loader.dataset)} ({100. * idx / len(train_loader):.3f}%)]\t'
                  f'Loss: {loss.item():.3f}\t\t', end='')
    print()


class NSTrainSet(Dataset):
    """
    Negative sampling completamente aleatorio.

    PARCHADO (TFM-NetActivity, graph_v3): el original pre-materializaba
    todo el dataset tileado (l*(1+neg) filas) en memoria de una sola vez
    via np.tile. Con l en el orden de cientos de millones esto satura
    la RAM antes de que arranque el entrenamiento (OOM). Esta version
    guarda solo las l muestras positivas (sample) y genera el nodo
    negativo al vuelo en __getitem__, sin cambiar la distribucion de
    muestreo ni la interfaz que espera train().
    """

    def __init__(self, sample, node_size, neg=5):
        """
        :param node_size: numero de nodos
        :param neg: numero de negativos por muestra positiva
        :param sample: retorno de HIN.sample(), (start_node, end_node, path_id)
        """
        print('init training dataset...')
        self.sample = np.asarray(sample, dtype=np.int64)  # shape (l, 3)
        self.l = len(self.sample)
        self.node_size = node_size
        self.neg = neg
        self.length = self.l * (1 + neg)
        print('finished')

    def __getitem__(self, index):
        base_idx = index % self.l
        start, end, path = self.sample[base_idx]
        if index < self.l:
            y = 1.0
        else:
            y = 0.0
            end = np.random.randint(0, self.node_size - 1)
        x = torch.LongTensor([start, end, path])
        y = torch.tensor(y, dtype=torch.float32)
        return x, y

    def __len__(self):
        return self.length
