import seaborn as sns
import pandas as pd
import numpy as np


def plot(df: pd.DataFrame, id_vars: list, value_vars: list):
    """
    Plots a data frame.

    :param df: dataframe to plot
    :type df: pd.DataFrame
    :param id_vars: list of column names for the x-axis
    :type id_vars: list or str
    :param value_vars: list of column names for the y-axis
    :type value_vars: list
    :return: melted data frame
    :rtype: pd.DataFrame
    """
    if not isinstance(id_vars, list):
        id_vars = [id_vars]
    data = pd.melt(df, id_vars=id_vars, value_vars=value_vars)

    sns.lineplot(x=id_vars[0], y="value",
                 hue="variable",
                 data=data)
    return data


def gradient(df: pd.DataFrame, xlabel: str, n=1):
    """
    Computes the gradient from a DataFrame

    :param df: dataframe to plot
    :type df: pd.DataFrame
    :param xlabel: column name for the x-axis
    :type xlabel: str
    :param n: the n-th derivative (default: 1)
    :type n: int
    :return: gradient as a DataFrame
    :rtype: pd.DataFrame
    """
    x = np.array(df[xlabel]).reshape(len(df[xlabel]), 1)

    first = np.gradient(df)[0]  # first derivative
    nth = first  # n-th derivative
    if n > 1:
        for _ in range(n - 1):
            nth = np.gradient(nth)[0]

    i = list(df.columns).index(xlabel)

    dx = first[:, i].reshape(len(x), 1)

    c = list(range(len(df.columns)))
    c.remove(i)
    dy = nth[:, c]

    dydx = dy / dx

    grad_df = pd.DataFrame(np.concatenate((dydx, x), axis=1), columns=df.columns)
    return grad_df
