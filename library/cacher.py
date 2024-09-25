import hashlib
import os
import numpy as np
import csv

from scipy.interpolate import griddata, interp1d
from scipy.optimize import curve_fit

class CachedResults:

    def __init__(self, cache_dir='cache'):
        self.cache_dir = cache_dir

    def create_cache_key(self, melt_comp, volatile_comp):
        melt_comp_tuple = tuple(sorted(melt_comp.items()))
        volatile_comp_tuple = tuple(sorted(volatile_comp.items()))
        key_string = f"{melt_comp_tuple}{volatile_comp_tuple}"
        return hashlib.md5(key_string.encode()).hexdigest()

    def create_cache_filename(self, melt_comp, volatile_comp):
        melt_comp_str = "_".join([f"{k}-{round(v, 2)}" for k, v in sorted(melt_comp.items())])
        volatile_comp_str = "_".join([f"{k}-{round(v, 2)}" for k, v in sorted(volatile_comp.items())])
        return f"melt_{melt_comp_str}__volatile_{volatile_comp_str}.csv"

    def save_to_cache(self, melt_comp, volatile_comp, T, P, fO2, O_abun, mb_eq):
        os.makedirs(self.cache_dir, exist_ok=True)
        cache_filename = self.create_cache_filename(melt_comp, volatile_comp)
        cache_path = os.path.join(self.cache_dir, cache_filename)

        # Retrieve all existing results
        all_results = self.retrieve_all_cached_results(melt_comp, volatile_comp)
        # Add the new result
        all_results.append((T, P, fO2, O_abun, mb_eq))        

        # Sort results by T first and then by P
        all_results.sort(key=lambda x: (x[0], x[1]))

        # Write sorted results to the cache file
        with open(cache_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['T', 'P', 'fO2', 'O_abun', 'mb_eq'])
            writer.writerows(all_results)

    def retrieve_all_cached_results(self, melt_comp, volatile_comp):
        cache_filename = self.create_cache_filename(melt_comp, volatile_comp)
        cache_path = os.path.join(self.cache_dir, cache_filename)
        
        if os.path.exists(cache_path):
            with open(cache_path, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                results = [(float(row['T']), float(row['P']), float(row['fO2']),\
                            float(row['O_abun']), float(row['mb_eq'])) for row in reader]
            return results
        return []

    def get_cached_value(self, cached_results, T, P):
        for cached_T, cached_P, cached_fO2, cached_O_abun, cached_mb_eq in cached_results:
            if cached_T == T and cached_P == P:
                return cached_fO2, cached_O_abun, cached_mb_eq
        return None

    def linear_extrapolation_1d(self, x, y, target):
        def linear_model(x, a, b):
            return a * x + b
        
        popt, _ = curve_fit(linear_model, x, y)
        return linear_model(target, *popt)

    def interpolate_or_extrapolate_results(self, results, T, P):
        if len(results) < 4:
            return None
        
        points = np.array([(result[0], result[1]) for result in results])
        values = np.array([result[2] for result in results])
        
        T_unique = np.unique(points[:, 0])
        print('T_unique',T_unique)
        print(len(T_unique))
        P_unique = np.unique(points[:, 1])
        
        T_min, T_max = points[:, 0].min(), points[:, 0].max()
        P_min, P_max = points[:, 1].min(), points[:, 1].max()
        
        if T_min <= T <= T_max and P_min <= P <= P_max:
            if len(T_unique) > 1 and len(P_unique) > 1:
                print('Point falls between existing cache values.') 
                print('Using interpolation for first estimate of fO2.')
                estimated_value = griddata(points, values, (T, P), method='linear')
                print(f'Estimated value: {estimated_value}')
                return estimated_value
            elif len(T_unique) > 1:
                print('Interpolating along P dimension.')
                f = interp1d(points[:, 0], values, kind='linear', fill_value="extrapolate")
                return f(T)
            elif len(P_unique) > 1:
                print('Interpolating along T dimension.')
                f = interp1d(points[:, 1], values, kind='linear', fill_value="extrapolate")
                return f(P)
        
        elif (T_min - 200 <= T <= T_max + 200) and (10**(np.log10(P_min) - 1) <= P <= 10**(np.log10(P_max) + 1)):
            print('Point falls within 200 K and 1 order of magnitude of cached values.')
            print('Using extrapolation for first estimate of fO2.')
            if len(T_unique) > 1 and len(P_unique) > 1:
                interpolator = LinearNDInterpolator(points, values)
                estimated_value = interpolator(T, P)
                print(f'Estimated value: {estimated_value}')
                return estimated_value
            elif len(P_unique) > 1:
                print('Extrapolating along P dimension.')
                return self.linear_extrapolation_1d(points[:, 0], values, T)
            elif len(T_unique) > 1:
                print('Extrapolating along T dimension.')
                return self.linear_extrapolation_1d(points[:, 1], values, P)
        
        return None