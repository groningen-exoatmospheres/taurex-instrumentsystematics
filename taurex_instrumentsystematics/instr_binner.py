import numpy as np
from taurex.binning.binner import Binner
from taurex.util.util import compute_bin_edges, create_grid_res, wnwidth_to_wlwidth
from taurex import OutputSize
from taurex.binning import FluxBinner
from astropy.io import fits

class FluxBinnerConv(Binner):
    """
    Bins to a wavenumber grid given by ``wngrid`` using a
    more accurate method that takes into account the amount
    of contribution from each native bin. This method also
    handles cases where bins are not continuous and/or
    overlapping.

    Parameters
    ----------

    wngrid: :obj:`array`
        Wavenumber grid

    wngrid_width: :obj:`array`, optional
        Must have same shape as ``wngrid``
        Full bin widths for each wavenumber grid point
        given in ``wngrid``. If not provided then
        this is automatically computed from ``wngrid``.

    """
    def __init__(self, wlgrids, wlgrid_widths=None, 
                 broadening_profiles = None,
                 profile_type = 'stsci_fits',
                 max_wlbroadening = None,
                 factor_cut = 5,
                 wlres = 15000):
        super().__init__()
        
        self._wlgrids = wlgrids
        self._wlgrid_widths = wlgrid_widths
        self._broadening_profiles = broadening_profiles
        self._profile_type = profile_type
        self._max_wlbroadening = max_wlbroadening
        self._factor_cut = factor_cut
        self._wlres = wlres
        
        self._wlgrid = np.concatenate(self._wlgrids)
        self._wlgrid_width = np.concatenate(self._wlgrid_widths)
        
        if self._wlgrid_widths is None:
            raise WlGridWidthsNotProvidedError('Use compute_bin_edges function to provide wngrid widths')
            
        self.binners = []
        self.sorters = []
        for i, w in enumerate(self._wlgrids):
            s = w.argsort()
            self.sorters.append(s)
            self.binners.append(FluxBinner(w[s], self._wlgrid_widths[i][s]))

        self._wlgrid = np.concatenate(self._wlgrids)
        self._wlgrid_width = np.concatenate(self._wlgrid_widths)
        ###sort_grid = wngrid.argsort()
        ###self._wngrid = wngrid[sort_grid]
        ###self._wngrid_width = wngrid_width
###
        ###if self._wngrid_width is None:
        ###    self._wngrid_width = compute_bin_edges(self._wngrid)[-1]
        ###elif hasattr(self._wngrid_width, '__len__'):
        ###    if len(self._wngrid_width) != len(self._wngrid):
        ###        raise ValueError('Wavenumber width should be signel value or '
        ###                         'same shape as wavenumber grid')
        ###    self._wngrid_width = wngrid_width[sort_grid]
###
        ###if not hasattr(self._wngrid_width, '__len__'):
        ###    self._wngrid_width = np.ones_like(self._wngrid)*self._wngrid_width
        
        if self._profile_type == 'stsci_fits':
            self._profiles, self._grid_fbs = self.load_stsci_fits(self._broadening_profiles)
        
    def load_stsci_fits(self, files):
        profiles = []
        grid_fbs = []
        for i, f in enumerate(files):
            try:
                with fits.open(f) as hdu:
                    science_data = hdu[1].data
                wl = science_data['WAVELENGTH'][:]
                R = science_data['R'][:]
            except:
                science_data = np.loadtxt(f)
                wl = science_data[:,0]
                R = science_data[:,1]
            std = wl/R/2 ### /2 to transform widths in std
            new_grid = create_grid_res(self._wlres, self._wlgrids[i][0]-10*std[0], self._wlgrids[i][-1]+10*std[-1])
            grid_fbs.append(FluxBinner(new_grid[:,0], new_grid[:,1]))
            s = np.interp(new_grid[:,0], wl, std, left=std[0], right=std[-1], period=None)
            profiles.append(np.clip(s, a_min=1e-20, a_max=self._max_wlbroadening))
        return profiles, grid_fbs
    
    def gaussian(self, x, mean, std):
        return 1.0 / (np.sqrt(2.0 * np.pi) * std) * np.exp(-np.power((x - mean) / std, 2.0) / 2)
            
    def low_res_convolved(self, Y, prof):
        #native_grid, final_flux, ers, widths = O
        #Y = fb.bindown(native_grid, final_flux, error=ers, grid_width=widths)
        #X = 10000/Y[0], Y[1]

        X_conv = np.zeros(len(Y[0]))
        for i in range(len(Y[0])):
            #std = np.clip(self._wlbroadening_default + self._wlbroadening_width*X[0][i] + self._wlbroadening_width2*X[0][i]**2, a_min=1e-20, a_max=self._max_wlbroadening)
            std = prof[i]
            if i != len(Y[0])-1:
                windows = np.maximum(1, int(self._factor_cut*std/np.abs(Y[0][i+1]-Y[0][i])))
            else:
                windows = np.maximum(1, int(self._factor_cut*std/np.abs(Y[0][i]-Y[0][i-1])))
            imin = np.maximum(0, i-windows)
            imax = np.minimum(len(Y[0]), i+windows+1)
            x = Y[0][imin:imax]
            y = self.gaussian(x, Y[0][i], std)
            X_conv[imin:imax] += Y[1][i]*y/np.sum(y)

        return Y[0], X_conv, Y[2], Y[3]

    def bindown(self, wngrid, spectrum, grid_width=None, error=None, in_wavenumber = True):
        if isinstance(spectrum, np.ndarray):
            return self.bindown2d(wngrid, spectrum, grid_width=grid_width, error=error, in_wavenumber = in_wavenumber)
        else:
            return self.bindown1d(wngrid, spectrum, grid_width=grid_width, error=error, in_wavenumber = in_wavenumber)
        #elif spectrum.ndim == 2:
            

    def bindown1d(self, wngrid, spectrum, grid_width=None, error=None, in_wavenumber = True):
        if in_wavenumber:
            if grid_width is not None:
                grid_width = wnwidth_to_wlwidth(wngrid, grid_width)[::-1]
            if error is not None:
                error = error[::-1]

            O_master = (10000/wngrid[::-1], spectrum[::-1], error, grid_width) ## errror and grid widths are already reversed and in correct units.
        else:
            O_master = (wngrid, spectrum, error, grid_width)
        wls = []
        sps = []
        ers = []
        wws = []
        for i, b in enumerate(self.binners):
            if self._profile_type == 'stsci_fits':
                o = self._grid_fbs[i].bindown(O_master[0], O_master[1], error=O_master[3], grid_width=O_master[2])
                #self.SAVE1 = o
                o = self.low_res_convolved(o, self._profiles[i])
                #self.SAVE2 = o
            else:
                o = O_master
            o = b.bindown(o[0], o[1], grid_width=o[3], error=o[2])
            #self.SAVE3 = o
            
            wls.append(o[0])
            sps.append(o[1])
            ers.append(o[2])
            wws.append(o[3])
        wlgrid = np.concatenate(wls)
        wlgrid_width = np.concatenate(wws)
        if error == None:
            res_ers = None
        else:
            res_ers = np.concatenate(ers)
   
        return wlgrid, np.concatenate(sps), res_ers, wlgrid_width

    def bindown2d(self, wngrid, spectrum, grid_width=None, error=None, in_wavenumber = True):
        Nspec = len(spectrum)

        for i, b in enumerate(self.binners):

            if in_wavenumber:
                if grid_width is not None:
                    grid_width = wnwidth_to_wlwidth(wngrid[i], grid_width[i])[::-1]
                if error is not None:
                    error = error[i][::-1]

                O_master = (10000/wngrid[i][::-1], spectrum[i][::-1], error, grid_width) ## errror and grid widths are already reversed and in correct units.
            else:
                if grid_width is not None:
                    grid_width = grid_width[i]
                if error is not None:
                    error = error[i]
                O_master = (wngrid[i], spectrum[i], error, grid_width)
            wls = []
            sps = []
            ers = []
            wws = []
        
            if self._profile_type == 'stsci_fits':
                o = self._grid_fbs[i].bindown(O_master[0], O_master[1], error=O_master[3], grid_width=O_master[2])
                #self.SAVE1 = o
                o = self.low_res_convolved(o, self._profiles[i])
                #self.SAVE2 = o
            else:
                o = O_master
            o = b.bindown(o[0], o[1], grid_width=o[3], error=o[2])
            #self.SAVE3 = o
            
            wls.append(o[0])
            sps.append(o[1])
            ers.append(o[2])
            wws.append(o[3])
        wlgrid = np.concatenate(wls)
        wlgrid_width = np.concatenate(wws)
        if error == None:
            res_ers = None
        else:
            res_ers = np.concatenate(ers)
   
        return wlgrid, np.concatenate(sps), res_ers, wlgrid_width

    def generate_spectrum_output(self, model_output,
                                 output_size=OutputSize.heavy):

        output = {}
        wngrid, flux, tau, extra = model_output
        output['native_wngrid'] = wngrid
        output['native_wlgrid'] = 10000/wngrid
        output['native_spectrum'] = flux
        output['native_wnwidth'] = compute_bin_edges(wngrid)[-1]

        output['binned_wngrid'] = 10000/self._wlgrid
        output['binned_wlgrid'] = self._wlgrid
        #output['binned_wnwidth'] = self._wngrid_width
        output['binned_wlwidth'] = self._wlgrid_width
        output['binned_spectrum'] = self.bindown(wngrid, flux)[1]
        if output_size > OutputSize.lighter:
            #output['binned_tau'] = self.bindown(wngrid, tau)[1]
            if output_size > OutputSize.light:
                output['native_tau'] = tau

        return output
