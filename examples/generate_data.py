from tf_pwa.config_loader import ConfigLoader
from tf_pwa.data import data_to_numpy, flatten_dict_data
import numpy as np
from tf_pwa.config import set_config

def main():
    set_config("polar", False)
    config = ConfigLoader("config.yml")
    config.set_params("toy_params.json")
    print(config.get_params())
    phsp = config.generate_phsp(10000)
    np.savez("phsp.npz", **data_to_numpy(flatten_dict_data(phsp)))
    toy = config.generate_toy(1000)
    np.savez("toy.npz", **data_to_numpy(flatten_dict_data(toy)))
    amp = config.get_amplitude()
    print(amp(toy)[:10])

if __name__ == "__main__":
    main()
