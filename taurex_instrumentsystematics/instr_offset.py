import numpy as np
import copy
from taurex.data.spectrum import ArraySpectrum, BaseSpectrum
from taurex.util.util import wnwidth_to_wlwidth, compute_bin_edges, create_grid_res
from taurex.log import Logger
from taurex.data.fittable import fitparam
from .instr_binner import FluxBinnerConv


class OffsetSpectraCont(BaseSpectrum):
    def __init__(self, path_spectra = [], offsets = [], slopes = [], slope_type ='linear',
                    broadening_profiles = [], wlshift = 0.0, max_wlbroadening = 0.1, factor_cut = 5, wlres = 15000,):
        
        super().__init__(self.__class__.__name__)
        
        self._wlshift = wlshift
        self._broadening_profiles = broadening_profiles
        self._max_wlbroadening = max_wlbroadening
        self._factor_cut = factor_cut
        self._wlres = wlres
        
        #Logger.__init__(self,'MultiSpectra')
        self.path_spectra = path_spectra
        self.slope_type = slope_type
        self.offsets = offsets
        self.slopes = slopes
        if len(self.offsets) == 0:
            #self.info('offsets are not init, set to 0')
            self.offsets = [0.0]*len(self.path_spectra)
        elif len(self.offsets) != len(self.path_spectra):
            #self.error('offsets and spectra are of different dimension')
            raise NotImplementedError('Check dimensions of the parameters path_spectra and offsets')
        if len(self.slopes) != len(self.path_spectra): 
            self.slopes = [0.0]*len(self.path_spectra)
        
        self._raw = []
        for i, s in enumerate(self.path_spectra):
            data = np.loadtxt(s)
            #data = data[data[:, 0].argsort(axis=0)[::-1]]
            if data.ndim == 1:
                self._raw.append(np.array([data]))
                pass
            elif data.ndim ==2:
                self._raw.append(data)
            else:
                raise NotImplementedError('Spectra are not in a 2D array format...')
        
        #self.spectrum
        self._raw = self._sort_spectra(self._raw)
        
        spec = copy.deepcopy(self._raw)
        for i, o in enumerate(self.offsets):
            spec[i][:,1] = spec[i][:,1]+o
            if i == 0:
                self._obs_spectrum = spec[i]
            else:
                self._obs_spectrum = np.vstack((self._obs_spectrum,spec[i]))
        #self.sort = self._obs_spectrum[:, 0].argsort(axis=0)[::-1]
        
        #super().__init__(self._obs_spectrum)
        
        #self._obs_spectrum = spectrum
        self._bin_widths = None
        self._bin_edges = None

        #self._sort_spectrum()
        self._process_spectrum()

        #self._wnwidths = wnwidth_to_wlwidth(self.wavelengthGrid,
        #                                    self._bin_widths)
        
        self.generate_offset_fitting_params()

    #def _sort_spectrum(self):
    #    self._obs_spectrum = \
    #        self._obs_spectrum[self._obs_spectrum[:, 0].argsort(axis=0)[::-1]]
    
    def _sort_spectra(self, spec):
        final_spec = []
        for i, s in enumerate(spec):
            idxs = s[:,0].argsort()
            X = s[idxs,:]
            final_spec.append(X)
        return final_spec

    def _process_spectrum(self):
        """
        Seperates out the observed data, error, grid and binwidths
        from array. If bin widths are not present then they are
        calculated here
        """
        self._raw_bin_edges = []
        self._raw_bin_widths = []
        self._raw_wnwidths = []
        for i,s in enumerate(self._raw):
            if s.shape[1] == 3:
                bin_edges, bin_widths = compute_bin_edges(s[:,0])
                self._raw_bin_edges.append(bin_edges)
                self._raw_bin_widths.append(bin_widths)
            elif s.shape[1] == 4:
                bin_widths = s[:,3]
                self._raw_bin_widths.append(bin_widths)
                obs_wl = s[:,0][::-1]
                obs_bw = bin_widths[::-1]

                bin_edges = np.zeros(shape=(len(bin_widths)*2,))

                bin_edges[0::2] = obs_wl - obs_bw/2
                bin_edges[1::2] = obs_wl + obs_bw/2

                self._raw_bin_edges.append(bin_edges[:])
            else:
                raise ErrorOfObsSpectrumShape
            
            self._raw_wnwidths.append(wnwidth_to_wlwidth(s[:,0],
                                            bin_widths))
        self._wnwidths = np.concatenate(self._raw_wnwidths)
        self._bin_edges = np.concatenate(self._raw_bin_edges)

    @property
    def rawData(self):
        """Data read from file"""
        return self._obs_spectrum

    #@property
    #def spectrum(self):
    #    """The spectrum itself"""
    #    return self._obs_spectrum[:, 1]

    @property
    def wavelengthGrid(self):
        """Wavelength grid in microns"""
        return self.rawData[:, 0]

    @property
    def wavenumberGrid(self):
        """Wavenumber grid in cm-1"""
        return 10000/self.wavelengthGrid

    @property
    def binEdges(self):
        """ Bin edges"""
        return 10000/self._bin_edges[:]

    @property
    def binWidths(self):
        """bin widths"""
        return self._wnwidths[:]

    @property
    def errorBar(self):
        """ Error bars for the spectrum"""
        return self.rawData[:, 2]
        
    def create_binner(self):
        """
        Creates the appropriate binning object
        """
        #from taurex.binning import FluxBinner
        self._raw_wngrid = [10000/r[:,0] for r in self._raw]
        self._raw_wlgrid = [r[:,0] for r in self._raw]

        return FluxBinnerConv(wlgrids=self._raw_wlgrid,
                            wlgrid_widths=self._raw_bin_widths,
                            broadening_profiles = self._broadening_profiles,
                            max_wlbroadening = self._max_wlbroadening,
                            factor_cut = self._factor_cut,
                            wlres = self._wlres)
    
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
        spec = copy.deepcopy(self._raw)
        for i, o in enumerate(self.offsets):
            if (self.slope_type == 'linear') or (self.slope_type is None):
                spec[i][:,1] = spec[i][:,1] + o + self.slopes[i]* (spec[i][:,0]- np.mean(spec[i][:,0]))
            elif self.slope_type == 'log':
                spec[i][:,1] = spec[i][:,1] + o + self.slopes[i]* 10**(spec[i][:,0]- np.mean(spec[i][:,0]))
            if i == 0:
                X = spec[i][:,1]
            else:
                X = np.concatenate( (X,spec[i][:,1]))
        return X[:]
        #X = X[self.sort]
        #self._obs_spectrum[:,1] = X
        #return self._obs_spectrum[:,1]
    
    @classmethod
    def input_keywords(self):
        return ['spectra_instr', 'observation_instr']