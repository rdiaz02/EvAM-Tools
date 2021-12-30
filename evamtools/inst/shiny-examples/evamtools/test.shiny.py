import pdb
import unittest
from selenium import webdriver
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys
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
    
    def _get_table_info(self):
        row_table = self.driver.find_elements_by_css_selector("#csd_freqs tbody>tr")
        row_text = [i.text for i in row_table]
        return row_text

    def _modify_genotype(self, genotype = None, freq = None):
        if genotype:
            for gene in genotype.split(", "):
                input_gene = self.driver.find_element_by_css_selector(
                    f"#genotype input[value={gene}]")
                input_gene.find_element_by_xpath('..').click()
        
        if freq != None:
            genotype_freq = self.driver.find_element_by_css_selector("#genotype_freq")
            genotype_freq.clear()
            genotype_freq.send_keys(freq)

        self.driver.find_element_by_css_selector("#add_genotype").click()
        sleep(1)

    def _select_tab(self, input_type, name_dataset):
        csd_tab = self.driver.find_element_by_css_selector(f"#input2build .radio input[value={input_type}]")
        csd_tab.click()
        dataset_tab = self.driver.find_element_by_css_selector(f"#select_csd .radio input[value={name_dataset}]")
        dataset_tab.click()
        sleep(1)

    def _process_genotype_table(self, data):
        genotypes = {}

        # pdb.set_trace()
        for i in data:
            key, val = i.replace(", ", "").split(" ")
            genotypes[key] = int(val)
        
        return genotypes

    def _get_genotypes_dict(self):
        return self._process_genotype_table(self._get_table_info())


    ## TESTING BASIC FUNCIONALITY
    # def test_switch_tab(self):
        ## Go to CSD-LINEAR

        ## Go to DAG-AND

        ## Go to MATRIX-EXAMPLE

        ## Go back to CSD & check we are in LINEAR

        ## Go back to DAG & check we are in AND

        ## Go back to MATRIX & check we are in EXAMPLE

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
    def test_modiying_genotype_with_buttons(self):
        ## Initial status
        status = self._get_status()
        assert(status["number_of_genes"] == 3)
        assert(status["selected_dataset"] == "User")
        assert(status["selected_input2build"] == "csd")
        assert(status["gene_names"] == ["A", "B", "C"])

        table_info = self._get_table_info()
        assert(table_info[0] == "No data available in table")

        analysis_button = self.driver.find_element_by_css_selector("#analysis")

        ## Cannot add string frequencies
        self._modify_genotype("A", "asdf")
        sleep(0.5)
        table_info = self._get_table_info()
        assert(table_info[0] == "No data available in table")
        assert(analysis_button.is_enabled() == False)

        ## Cannot add empty frequencies
        self._modify_genotype("A")
        sleep(0.5)
        table_info = self._get_table_info()
        assert(table_info[0] == "No data available in table")
        assert(analysis_button.is_enabled() == False)

        ## Cannot add negative frequencies
        self._modify_genotype("A", -100)
        sleep(0.5)
        table_info = self._get_table_info()
        assert(table_info[0] == "No data available in table")
        assert(analysis_button.is_enabled() == False)

        ## Can add positive freq
        self._modify_genotype("A", 100)
        sleep(1)
        table_info = self._get_table_info()
        assert(table_info[0] == "A 100")
        assert(analysis_button.is_enabled())

        ## Freq changes based on available data 
        input_gene = self.driver.find_element_by_css_selector(
                    "#genotype input[value=A]")
        input_gene.find_element_by_xpath('..').click()
        current_freq = self.driver.find_element_by_css_selector("#genotype_freq").get_attribute("value")
        print(current_freq)
        assert(int(current_freq) == 100)

        ## Freq changes based on available data 
        input_gene = self.driver.find_element_by_css_selector(
                    "#genotype input[value=A]")
        input_gene.find_element_by_xpath('..').click()
        sleep(1)
        current_freq = self.driver.find_element_by_css_selector("#genotype_freq").get_attribute("value")
        print(current_freq)
        assert(current_freq == '')

        ## Remove genotype when adding a 0 freq
        sleep(0.5)
        self._modify_genotype(genotype = "A", freq = 0)
        sleep(1)
        table_info = self._get_table_info()
        assert(table_info[0] == "No data available in table")

        ## Not specifying genes add WT freq
        self._modify_genotype(freq = 100)
        sleep(1)
        table_info = self._get_table_info()
        assert(table_info[0] == "WT 100")
        assert(analysis_button.is_enabled())

    def test_modify_table(self):
        ## Select AND dataset
        self._select_tab("csd", "AND")

        ## Check Status
        status = self._get_status()
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "AND")
        assert(status["selected_input2build"] == "csd")
        assert(status["gene_names"] == ["A", "B", "C", "D"])

        ## Table before modification
        prev_table_info = self._get_table_info()
        prev_genotypes = self._process_genotype_table(prev_table_info)

        ## Double click on first freq
        first_freq = self.driver.find_elements_by_css_selector("#csd_freqs tbody>tr>td")[1]
        actions = ActionChains(self.driver)
        actions.double_click(first_freq).perform()
        # pdb.set_trace()
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.RIGHT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).perform()

        ## Set freqs (remove some genotypes)
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.LEFT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys(200).perform()
        actions.send_keys(Keys.TAB).perform()
        for i in [50, 150, 0]:
            actions.send_keys(i)
            actions.send_keys(Keys.TAB).perform()

        ## Save changes
        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()

        ## Check Status
        sleep(1)
        out_table_info = self._get_table_info()
        out_genotypes = self._process_genotype_table(out_table_info)

        assert(out_genotypes["WT"] == 200)
        assert(out_genotypes["A"] == 50)
        assert(out_genotypes["AB"] == 150)
        assert(("ABC" in out_genotypes.keys()) == False)
        assert(out_genotypes["AC"] == prev_genotypes["AC"])
        assert(out_genotypes["ABCD"] == prev_genotypes["ABCD"])

    def test_change_gene_number_CSD(self):
        ## Select Linear & save data
        self._select_tab("csd", "Linear")
        initial_genotypes = self._get_table_info()
        initial_genotypes_dict = self._process_genotype_table(initial_genotypes)

        ## Change to 5 genes
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()

        sleep(1)

        ## Add new genotype & save data
        status = self._get_status()
        assert(status["gene_names"] == ["A", "B", "C", "D", "E", "F", "G"])
        self._modify_genotype("A, B, C, D, E", 100)
        new_genotypes_dict = self._get_genotypes_dict()
        assert(new_genotypes_dict["ABCDE"] == 100)

        ## Set to 4 genes & check we lose some data
        move.click_and_hold(slider_input).move_by_offset(-200, 0).release().perform()
        sleep(1)
        status = self._get_status()
        assert(status["gene_names"] == ["A", "B", "C", "D"])
        new_genotypes_dict = self._get_genotypes_dict()
        assert("ABCDE" not in new_genotypes_dict.keys())        

        ## Set to 5 genes & check we recover some data
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()
        sleep(1)
        new_genotypes_dict = self._get_genotypes_dict()
        assert(new_genotypes_dict["ABCDE"]== 100)

    def test_change_gene_names_CSD(self):
        ## Select Linear & save data
        self._select_tab("csd", "Linear")
        initial_genotypes = self._get_table_info()
        initial_genotypes_dict = self._process_genotype_table(initial_genotypes)

        base_names = ["A", "B", "C", "D"]
        new_names = ["A1", "B2", "C3", "D4"]

        expected_genotypes_dict = {}
        for k, v in initial_genotypes_dict.items():
            for (base_gene_name, new_gene_name) in zip(base_names, new_names):
                k = k.replace(base_gene_name, new_gene_name)
            expected_genotypes_dict[k] = v

        ## Change gene names
        self.driver.find_element_by_css_selector("#change_gene_names").click()
        input_new_gene_names = self.driver.find_element_by_css_selector("#new_gene_names")
        input_new_gene_names.clear()
        input_new_gene_names.send_keys("A1, B2, C3, D4")
        self.driver.find_element_by_css_selector("#action_gene_names").click()

        actions = ActionChains(self.driver)
        actions.move_to_element_with_offset(input_new_gene_names, -500, 500)
        actions.click()
        actions.perform()
        sleep(1)

        ## Get gene names
        status = self._get_status()
        assert(status["gene_names"], new_names)

        ## Get genotypes freqs
        genotype_info = self._get_table_info()
        genotype_info_dict = self._process_genotype_table(genotype_info)

        assert(genotype_info_dict == expected_genotypes_dict)

    def test_save_data_set_CSD(self):
        ## Selecting LINEAR
        self._select_tab("csd", "Linear")
        initial_genotypes = self._get_table_info()
        initial_genotypes_dict = self._process_genotype_table(initial_genotypes)

        ## Introducing modifications
        changes = {"A, B": 1000, "A": 500}

        for (genotype, freq) in changes.items():
            self._modify_genotype(genotype, freq)
            sleep(1.5)

        ## Saving with new name
        new_dataset_name = "SELENIUM_ft"

        save_button = self.driver.find_element_by_id("save_csd_data")
        assert(save_button.is_enabled() == False)
        dataset_name = self.driver.find_element_by_css_selector("input#dataset_name")
        dataset_name.clear()
        dataset_name.send_keys(new_dataset_name)
        sleep(0.5)
        assert(save_button.is_enabled() == True)
        save_button.click()

        sleep(0.5)
        ## Check status
        status = self._get_status()
        assert(status["selected_dataset"] == new_dataset_name)
        assert(status["selected_input2build"] == "csd")

        new_genotypes = self._get_table_info()
        new_genotypes_dict = self._process_genotype_table(new_genotypes)

        for genotype in new_genotypes_dict.keys():
            changed_genotypes = {}
             
            for k, v in changes.items():
                changed_genotypes[k.replace(", ", "")] = v

            if (genotype in changed_genotypes.keys()):
                assert(new_genotypes_dict[genotype] != initial_genotypes_dict[genotype])
                assert(new_genotypes_dict[genotype] == changed_genotypes[genotype])
            else:
                assert(new_genotypes_dict[genotype] == initial_genotypes_dict[genotype])

        ## Check dataset restoration
        self._select_tab("csd", "Linear")
        final_genotypes = self._get_table_info()

        assert(final_genotypes == initial_genotypes)

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
