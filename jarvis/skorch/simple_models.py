"""
Skorch with JARVIS-ML

From https://github.com/skorch-dev/skorch/blob/efeead6f1c361a4c12d35e78dd84eb7fe1c5232e/notebooks/Basic_Usage.ipynb
"""
from skorch import NeuralNetRegressor
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline


class RegressorModule(nn.Module):
    def __init__(self, num_units=10, nonlin=F.relu):
        super(RegressorModule, self).__init__()
        self.num_units = num_units
        self.nonlin = nonlin

        self.dense0 = nn.Linear(20, num_units)
        self.nonlin = nonlin
        self.dense1 = nn.Linear(num_units, 10)
        self.output = nn.Linear(10, 1)

    def forward(self, X, **kwargs):
        X = self.nonlin(self.dense0(X))
        X = F.relu(self.dense1(X))
        X = self.output(X)
        return X


def simple_regression(X=[], Y=[], plot=False, preprocess=True):
    """
  Quickly train simple regression models without hyperparameter-optimization
  Args:
      X: input features
      Y: Target data
      plot: whether to make a parity plot with ML models
      preprocess: whether to apply standard preprocessing techniques
  Returns:
  """

    net_regr = NeuralNetRegressor(
        RegressorModule,
        max_epochs=20,
        lr=0.1,
        #     device='cuda',  # uncomment this to train with CUDA
    )
    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        X, Y, jid, random_state=1, test_size=0.1
    )
    pipe = Pipeline(
        [
            ("stdscal", StandardScaler()),
            ("vart", VarianceThreshold(1e-4)),
            ("est", net_regr),
        ]
    )
    if preprocess == True:
        model = pipe
    else:
        model = net_regr
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    reg_sc = regr_scores(y_test, pred)
    if plot == True:
        plt.plot(reg_sc["pred"], reg_sc["test"], ".", label=str(type(i).__name__)[0:4])
    print(type(i).__name__, round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))
    if plot == True:
        plt.legend()
        plt.xlabel("DFT")
        plt.ylabel("ML")


if __name__ == "__main__":
    X, Y, jid = jdata(prop="exfoliation_en")
    simple_regression(x=X, y=Y)
