import pdb
import unittest
from selenium import webdriver
from time import sleep

class evamtools(unittest.TestCase):
    def setUp(self):
        self.driver = webdriver.Chrome('chromedriver')
        self.driver.implicitly_wait(1)
        self.driver.maximize_window()
        self.driver.get("http://127.0.0.1:3000/")
    
    def _get_status(self, type_of_input = "csd"):
        number_of_genes = int(self.driver.find_element_by_css_selector("#gene_number").get_attribute("value"))
        selected_input2build = self.driver.find_element_by_css_selector("#input2build input:checked").get_attribute("value")
        selected_dataset = self.driver.find_element_by_css_selector("#select_csd input:checked").get_attribute("value")

        if(type_of_input == "csd"):
            gene_names = self.driver.find_elements_by_css_selector("#genotype .checkbox input[type=checkbox]")
            gene_names = [i.get_attribute("value") for i in gene_names]
        elif(type_of_input == "dag"): 
            gene_names = self.driver.find_elements_by_css_selector("#dag_from input[type=radio]")
            gene_names = [i.get_attribute("value") for i in gene_names][1:]
        elif(type_of_input == "matrix"):
            gene_names = self.driver.find_elements_by_css_selector("#thetas_table th")
            gene_names = [i.text for i in gene_names][1:]
        
        out = {
           "number_of_genes": number_of_genes,
           "selected_input2build": selected_input2build,
           "gene_names": gene_names,
           "selected_dataset": selected_dataset,
        }

        return out
    
    def _get_error_message(self):
        return self.driver.find_element_by_css_selector("#shiny-modal .modal-body div")
        
    ## TESTING BASIC FUNCIONALITY
    # def test_switch_tab(self):
    #     csd_tab = self.driver.find_element_by_css_selector("#input2build .radio input[value=csd]")
    #     csd_tab.click()
    #     dag_tab = self.driver.find_element_by_css_selector("#input2build .radio input[value=dag]")
    #     dag_tab.click()
    #     matrix_tab = self.driver.find_element_by_css_selector("#input2build .radio input[value=matrix]")
    #     matrix_tab.click()

    def test_load_csv_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys("/home/pablo/Downloads/sample_csd/good_csd.csv")
        sleep(2)
        status = self._get_status()
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "good")
        assert(status["selected_input2build"] == "csd")
        assert(status["gene_names"] == ["A1", "B2", "C3", "D4"])

    def test_load_corrupt_csv_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys("/home/pablo/Downloads/sample_csd/bad_csd.csv")
        error_message = self._get_error_message().text
        assert(error_message, "Your csv data can not be loaded. Make sure it only contains 0 and 1.")
 
    def test_load_CSD_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys("/home/pablo/Downloads/sample_csd/FREQ_csd.RDS")
        sleep(2)
        status = self._get_status()
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "load_test")
        assert(status["selected_input2build"] == "csd")
        assert(status["gene_names"] == ["A", "B", "C", "D"])

    def test_load_DAG_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys("/home/pablo/Downloads/sample_csd/DAG_csd.RDS")
        sleep(2)
        status = self._get_status("dag")
        print(status)
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "DAG")
        assert(status["selected_input2build"] == "dag")
        assert(status["gene_names"] == ["A", "B", "C", "D"])

    def test_load_MATRIX_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys("/home/pablo/Downloads/sample_csd/MHN_csd.RDS")
        sleep(2)
        status = self._get_status("matrix")
        print(status)
        assert(status["number_of_genes"] == 3)
        assert(status["selected_dataset"] == "MHN")
        assert(status["selected_input2build"] == "matrix")
        assert(status["gene_names"] == ["A", "B", "C"])
    
    # ## TESTING CSD 
    # def test_add_genotype(self):
    #     pass
    
    # def test_remove_genotype(self):
    #     pass

    # def test_modify_table(self):
    #     pass

    # def test_change_gene_number_CSD(self):
    #     pass

    # def test_save_data_set_CSD(self):
    #     pass

    # def test_dag_pipeline(self):
    #     # 1 loading data
    #     # 2 modify dag
    #     # 3 modify lambda
    #     # 4 modify relathionship
    #     # 5 change gene names
    #     # 5.1 sample
    #     # 6 switch tabs
    #     # 7 go back
    #     # 8 save
    #     # 9 download
    #     # 10 run
    #     pass

    # ## TESTING DAG
    # def test_add_node(self):
    #     pass
    
    # def test_remove_node(self):
    #     pass

    # def test_make_cycle(self):
    #     pass

    # def test_repeated_relathionship(self):
    #     pass

    # def test_repeated_relathionship(self):
    #     pass

    # def test_change_names_DAG(self):
    #     pass

    # def test_change_gene_number_DAG(self):
    #     pass

    # def test_save_data_set_DAG(self):
    #     pass

    # def test_dag_pipeline(self):
    #     # 1 loading data
    #     # 2 modify dag
    #     # 3 modify lambda
    #     # 4 modify relathionship
    #     # 5 change gene names
    #     # 5.1 sample
    #     # 6 switch tabs
    #     # 7 go back
    #     # 8 save
    #     # 9 download
    #     # 10 run
    #     pass

    # ## TESTING MATRIX
    # def test_change_theta(self):
    #     pass
    
    # def test_repeated_relathionship(self):
    #     pass

    # def test_repeated_relathionship(self):
    #     pass

    # def test_change_names_MATRIX(self):
    #     pass

    # def test_change_gene_number_MATRIX(self):
    #     pass

    # def test_save_data_set_MATRIX(self):
    #     pass

    # def test_MATRIX_pipeline(self):
    #     # 1 loading data
    #     # 2 modify dag
    #     # 3 modify lambda
    #     # 4 modify relathionship
    #     # 5 change gene names
    #     # 5.1 sample
    #     # 6 switch tabs
    #     # 7 go back
    #     # 8 save
    #     # 9 download
    #     # 10 run
    #     pass

    ## TESTING RESULTS

if __name__ == "__main__":
    unittest.main()
