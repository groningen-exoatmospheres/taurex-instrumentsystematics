from taurex.data.spectrum import ArraySpectrum
from taurex.log import Logger
import copy
from taurex.data.fittable import fitparam
import numpy as np

class OffsetSpectra(ArraySpectrum):
    def __init__(self, path_spectra = [], offsets = [], slopes = [],error_scale=[], slope_type ='linear'):
        #Logger.__init__(self,'MultiSpectra')
        self.path_spectra = path_spectra
        self.slope_type = slope_type
        self.offsets = offsets
        self.escale = error_scale
        self.slopes = slopes
        if len(self.offsets) == 0:
            #self.info('offsets are not init, set to 0')
            self.offsets = [0.0]*len(self.path_spectra)
        elif len(self.offsets) != len(self.path_spectra):
            #self.error('offsets and spectra are of different dimension')
            raise NotImplementedError('Check dimensions of the parameters path_spectra and offsets')
        if len(self.slopes) != len(self.path_spectra): 
            self.slopes = [0.0]*len(self.path_spectra)
        if len(self.escale) != len(self.path_spectra): 
            self.escale = [1.0]*len(self.path_spectra)
        
        self._obs_spectra = []
        for i, s in enumerate(self.path_spectra):
            data = np.loadtxt(s)
            #data = data[data[:, 0].argsort(axis=0)[::-1]]
            if data.ndim == 1:
                self._obs_spectra.append(np.array([data]))
                pass
            elif data.ndim ==2:
                self._obs_spectra.append(data)
            else:
                raise NotImplementedError('Spectra are not in a 2D array format...')
        
        #self.spectrum
        
        spec = copy.deepcopy(self._obs_spectra)
        for i, o in enumerate(self.offsets):
            spec[i][:,1] = spec[i][:,1]+o
            if i == 0:
                self._obs_spectrum = spec[i]
            else:
                self._obs_spectrum = np.vstack((self._obs_spectrum,spec[i]))
        self.sort = self._obs_spectrum[:, 0].argsort(axis=0)[::-1]
        
        super().__init__(self._obs_spectrum)
        self.generate_offset_fitting_params()
        
    
    def generate_offset_fitting_params(self):

        bounds = [-0.001, 0.001]
        for idx, val in enumerate(self.offsets):

            point_num = idx+1
            param_name = 'Offset_{}'.format(point_num)
            param_latex = '$Offset_{}$'.format(point_num)
            def read_point(self, idx=idx):
                return self.offsets[idx]
            def write_point(self, value, idx=idx):
                self.offsets[idx] = value
            fget_point = read_point
            fset_point = write_point
            self.debug('FGet_location %s', fget_point)
            default_fit = False
            self.add_fittable_param(param_name, param_latex, fget_point,
                                    fset_point, 'linear', default_fit, bounds)

            param_name_escale = 'EScale_{}'.format(point_num)
            param_latex_escale = '$Escale_{}$'.format(point_num)
            def read_point_escale(self, idx=idx):
                return self.escale[idx]
            def write_point_escale(self, value, idx=idx):
                self.escale[idx] = value
            default_fit = False
            self.add_fittable_param(param_name_escale, param_latex_escale, read_point_escale,
                                    write_point_escale, 'linear', default_fit, bounds)

            param_name_slope = 'Slope_{}'.format(point_num)
            param_latex_slope = '$Slope_{}$'.format(point_num)
            def read_point_slope(self, idx=idx):
                return self.slopes[idx]
            def write_point_slope(self, value, idx=idx):
                self.slopes[idx] = value
            default_fit = False
            self.add_fittable_param(param_name_slope, param_latex_slope, read_point_slope,
                                    write_point_slope, 'linear', default_fit, bounds)
    
    @property
    def spectrum(self):
        spec = copy.deepcopy(self._obs_spectra)
        for i, o in enumerate(self.offsets):
            if (self.slope_type == 'linear') or (self.slope_type is None):
                spec[i][:,1] = spec[i][:,1] + o + self.slopes[i]* (spec[i][:,0]- np.mean(spec[i][:,0]))
            elif self.slope_type == 'log':
                spec[i][:,1] = spec[i][:,1] + o + self.slopes[i]* 10**(spec[i][:,0]- np.mean(spec[i][:,0]))
            if i == 0:
                X = spec[i][:,1]
            else:
                X = np.concatenate( (X,spec[i][:,1]))
        X = X[self.sort]
        #self._obs_spectrum[:,1] = X
        return X
    
    @property
    def errorBar(self):
        """ Error bars for the spectrum"""
        spec = copy.deepcopy(self._obs_spectra)
        for i, o in enumerate(self.escale):
            spec[i][:,2] = o*spec[i][:,2]
            if i == 0:
                X = spec[i][:,2]
            else:
                X = np.concatenate( (X,spec[i][:,2]))
        X = X[self.sort]
        #self._obs_spectrum[:,2] = X
        return X
    
    @classmethod
    def input_keywords(self):
        return ['spectra_w_offsets', 'observation_w_offsets']