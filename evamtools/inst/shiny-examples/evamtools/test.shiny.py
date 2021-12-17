import unittest
from selenium import webdriver

class evamtools(unittest.TestCase):
    def setUp(self):
        self.driver = None
        self.driver.implicitly_wait(30)
        self.driver.maximize_window()
        self.driver.get(None)
    
    ## TESTING BASIC FUNCIONALITY
    def test_switch_tab(self):
        pass

    def test_load_CSD_dataset(self):
        pass

    def test_load_DAG_dataset(self):
        pass

    def test_load_MATRIX_dataset(self):
        pass
    
    ## TESTING CSD 
    def test_add_genotype(self):
        pass
    
    def test_remove_genotype(self):
        pass

    def test_modify_table(self):
        pass

    def test_change_gene_number_CSD(self):
        pass

    def test_save_data_set_CSD(self):
        pass

    def test_dag_pipeline(self):
        # 1 loading data
        # 2 modify dag
        # 3 modify lambda
        # 4 modify relathionship
        # 5 change gene names
        # 5.1 sample
        # 6 switch tabs
        # 7 go back
        # 8 save
        # 9 download
        # 10 run
        pass

    ## TESTING DAG
    def test_add_node(self):
        pass
    
    def test_remove_node(self):
        pass

    def test_make_cycle(self):
        pass

    def test_repeated_relathionship(self):
        pass

    def test_repeated_relathionship(self):
        pass

    def test_change_names_DAG(self):
        pass

    def test_change_gene_number_DAG(self):
        pass

    def test_save_data_set_DAG(self):
        pass

    def test_dag_pipeline(self):
        # 1 loading data
        # 2 modify dag
        # 3 modify lambda
        # 4 modify relathionship
        # 5 change gene names
        # 5.1 sample
        # 6 switch tabs
        # 7 go back
        # 8 save
        # 9 download
        # 10 run
        pass

    ## TESTING MATRIX
    def test_change_theta(self):
        pass
    
    def test_repeated_relathionship(self):
        pass

    def test_repeated_relathionship(self):
        pass

    def test_change_names_MATRIX(self):
        pass

    def test_change_gene_number_MATRIX(self):
        pass

    def test_save_data_set_MATRIX(self):
        pass

    def test_MATRIX_pipeline(self):
        # 1 loading data
        # 2 modify dag
        # 3 modify lambda
        # 4 modify relathionship
        # 5 change gene names
        # 5.1 sample
        # 6 switch tabs
        # 7 go back
        # 8 save
        # 9 download
        # 10 run
        pass

    ## TESTING RESULTS
