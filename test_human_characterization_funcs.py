import human_characterisation_functions as hcf
import pytest
import pyabf

# def test_get_access_resistance():
#     filename_vc = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP220426_trial/22427001.abf'
#     ch = [2]
#     Rs, Rin = hcf.get_access_resistance(filename_vc, ch) 
#     assert Rs, Rin == (9.350773576, 37.91587646)

# def test_get_hyperpolar_param():
#     filename_char = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP220426_trial/22427003.abf'
#     channels = [2,3]
#     charact_dict = hcf.load_traces(filename_char)
#     inj=[-300,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,
#         600,700,800,900,1000,1100,1200,1300,1400]
#     tau_all, capacitance_all, mc_all, V65_all = hcf.get_hyperpolar_param(charact_dict, channels, inj)
#     assert capacitance_all == pytest.approx([445.3729258, 11.66120996])


def test_from_axon_file():
    fn = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP230523/23523050.abf' 

    data = hcf.from_axon_file(filepath = fn, channel = 6, unit ='pA', sweeps_delete = list(range(1,30)),
    first_point = 0, last_point = 5500)

    trace = pyabf.ABF(fn)
    trace.setSweep(sweepNumber = 0, channel = 6)
    data_cut = trace.sweepY
    data_cut = data_cut[5500:]

    assert (data == data_cut).all() 



