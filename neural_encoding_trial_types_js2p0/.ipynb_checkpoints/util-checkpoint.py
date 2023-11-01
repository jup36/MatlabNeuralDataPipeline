import os.path
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sparse
import glob
from pathlib import Path


class MatfileIO:
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.blockNums_dir = glob.glob(os.path.join(self.data_dir, 'blockNums_*.mat'))
        self.trialType_dir = glob.glob(os.path.join(self.data_dir, 'trialType_*.mat'))
        self.unitTimeBCtx_dir = glob.glob(os.path.join(self.data_dir, 'unitTimeBCtx_*.mat'))
        self.unitTimeBStr_dir = glob.glob(os.path.join(self.data_dir, 'unitTimeBStr_*.mat'))
        self.unitTimeBCg_dir = glob.glob(os.path.join(self.data_dir, 'unitTimeBCg_*.mat'))

    @classmethod
    def get_mat(cls, mat_dir):
        mat = []
        if mat_dir:
            mat = sio.loadmat(Path(mat_dir))
        return mat

    def get_var_as_list(self, dir_name, name_in_matlab):
        """
        This function uses the class method 'get_mat' to grab Matlab data of type 'str' or 'int' formatted as cell
        :param dir_name:
        :param name_in_matlab:
        :return:
        """
        _list = []
        if dir_name:
            _mat = self.get_mat(dir_name[0])
            var = np.squeeze(_mat[name_in_matlab])
            for i in range(len(var)):
                if name_in_matlab == 'ss_blNumbs':
                    _list.append(var[i][0][0])
                elif name_in_matlab == 'ss_trialType':
                    _list.append(var[i][:][0])
        return _list

    def get_sparse_as_dense(self, dir_name, name_in_matlab):
        """
        This function uses the class method 'get_mat' to grab Matlab data of type 'matrix (e.g. neuron-by-time)'
        formatted as cell
        sparse.csr_matrix.todense(mat['ss_unitTimeBCtx'][0][0])
        :param dir_name:
        :param name_in_matlab:
        :return:
        """
        _list = []
        if dir_name:
            _mat = self.get_mat(dir_name[0])
            var = np.squeeze(_mat[name_in_matlab])
            for i in range(len(var)):
                _list.append(sparse.csr_matrix.todense(var[i]))
        return _list

    def get_blockNums(self):
        return self.get_var_as_list(self.blockNums_dir, 'ss_blNumbs')

    def get_trialType(self):
        return self.get_var_as_list(self.trialType_dir, 'ss_trialType')

    def get_unitTimeBCtx(self):
        return self.get_sparse_as_dense(self.unitTimeBCtx_dir, 'ss_unitTimeBCtx')

    def get_unitTimeBStr(self):
        return self.get_sparse_as_dense(self.unitTimeBStr_dir, 'ss_unitTimeBStr')

    def get_unitTimeBCg(self):
        return self.get_sparse_as_dense(self.unitTimeBCg_dir, 'ss_unitTimeBCg')

    def extract_dataframe(self):
        # to dataframe
        df = pd.DataFrame()

        _bn = self.get_blockNums()
        df.insert(0, "blockNum", _bn)

        _tt = self.get_trialType()
        if len(_tt) == len(_bn):
            df.insert(1, "trialType", _tt)

        if len(self.unitTimeBCtx_dir) > 0:
            _ctx = self.get_unitTimeBCtx()
            if len(_ctx) == len(_bn):
                df.insert(df.shape[1], "ctx", _ctx)

        if len(self.unitTimeBStr_dir) > 0:
            _str = self.get_unitTimeBStr()
            if len(_str) == len(_bn):
                df.insert(df.shape[1], "str", _str)

        if len(self.unitTimeBCg_dir) > 0:
            _cg = self.get_unitTimeBCg()
            if len(_cg) == len(_bn):
                df.insert(df.shape[1], "cg", _cg)

        return df


class Bunch(dict):  # Bunch is a dictionary that supports attribute-style access, a la JavaScript
    """A subclass of dictionary with an additional dot syntax."""
    def __init__(self, *args, **kwargs):
        super(Bunch, self).__init__(*args, **kwargs)
        self.__dict__ = self

    def copy(self):
        """Return a new Bunch instance which is a copy of the current Bunch instance."""
        return Bunch(super(Bunch, self).copy())











