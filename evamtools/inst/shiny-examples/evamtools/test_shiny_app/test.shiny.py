import pdb
import os
import unittest
from unittest.case import expectedFailure
from selenium import webdriver
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys
from time import sleep

class evamtools_basics(unittest.TestCase):
    def setUp(self):
        self.driver = webdriver.Chrome('chromedriver')
        self.driver.implicitly_wait(1)
        self.driver.get("http://127.0.0.1:3000/")
        self.driver.maximize_window()
        self._go_to("csd_builder")
    
    def _get_status(self, type_of_input = "csd"):
        number_of_genes = int(self.driver.find_element_by_css_selector("#gene_number").get_attribute("value"))
        sleep(0.1)
        selected_input2build = self.driver.find_element_by_css_selector("#input2build input:checked").get_attribute("value")
        sleep(0.1)
        selected_dataset = self.driver.find_element_by_css_selector("#select_csd input:checked").get_attribute("value")
        sleep(0.1)

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
        return self.driver.find_element_by_css_selector("#shiny-modal .modal-body div").text
    
    def _get_table_info(self):
        row_table = self.driver.find_elements_by_css_selector("#csd_freqs tbody>tr")
        row_text = [i.text for i in row_table]
        return row_text

    def _go_to(self, where):
        self.driver.find_element_by_css_selector(f".nav a[data-value={where}]").click()
        sleep(0.5)

    def _modify_genotype(self, genotype = None, freq = None):
        if genotype:
            for gene in genotype.split(", "):
                input_gene = self.driver.find_element_by_css_selector(
                    f"#genotype input[value={gene}]")
                input_gene.find_element_by_xpath('..').click()
        sleep(0.5)
        if freq != None:
            genotype_freq = self.driver.find_element_by_css_selector("#genotype_freq")
            genotype_freq.clear()
            genotype_freq.send_keys(freq)

        sleep(0.5)
        self.driver.find_element_by_css_selector("#add_genotype").click()
        sleep(1)

    def _select_tab(self, input_type, name_dataset = None):
        self._scroll2top()
        csd_tab = self.driver.find_element_by_css_selector(f"#input2build .radio input[value={input_type}]")
        csd_tab.find_element_by_xpath('..').click()
        # csd_tab.click()
        sleep(0.7)
        if(name_dataset):
            dataset_tab = self.driver.find_element_by_css_selector(f"#select_csd .radio input[value={name_dataset}]")
            dataset_tab.find_element_by_xpath('..').click()
        sleep(0.7)
    
    def _select_result(self, cpm_name):
        dataset_tab = self.driver.find_element_by_css_selector(f"#cpm_list .radio input[value={cpm_name}]")
        dataset_tab.find_element_by_xpath('..').click()
        sleep(0.5)

    def _process_genotype_table(self, data):
        genotypes = {}

        for i in data:
            key, val = i.replace(", ", "").split(" ")
            genotypes[key] = int(val)
        
        return genotypes

    def _get_genotypes_dict(self):
        return self._process_genotype_table(self._get_table_info())

    def _change_gene_names(self, new_names):
        

        self.driver.find_element_by_css_selector("#change_gene_names").click()
        input_new_gene_names = self.driver.find_element_by_css_selector("#new_gene_names")
        input_new_gene_names.clear()
        input_new_gene_names.send_keys(", ".join(new_names))
        self.driver.find_element_by_css_selector("#action_gene_names").click()
        
        self.driver.find_element_by_css_selector(".modal-open .modal-footer button").click()
        sleep(1.5)
    
    def _add_edge(self, from_gene = None, to_gene = None):
        if(from_gene):
            self.driver.find_element_by_css_selector(f"#dag_from input[value={from_gene}]").click()
            sleep(0.5)
        
        if(to_gene):
            self.driver.find_element_by_css_selector(f"#dag_to input[value={to_gene}]").click()
            sleep(0.5)
        
        self.driver.find_element_by_css_selector("#add_edge").click()
        sleep(0.5)

    def _get_dag_info(self):
        row_table = self.driver.find_elements_by_css_selector("#dag_table tbody>tr")
        row_text = [i.text.split(" ") for i in row_table]
        return row_text

    def _get_theta_table(self):
        row_table = self.driver.find_elements_by_css_selector("#thetas_table tr")
        row_text = [i.text.split(" ") for i in row_table]
        return row_text

    def _scroll2top(self):
        self.driver.execute_script("window.scrollBy(0,-document.body.scrollHeight)")

## TESTING BASIC FUNCTIONALITY
class test_basic_functionality(evamtools_basics):
    # Testing basic input navigation
    def test_switch_tab(self):
        ## Go to CSD-LINEAR
        self._select_tab("csd", "Linear")

        ## Go to DAG-AND
        self._select_tab("dag", "AND")

        ## Go to MATRIX-EXAMPLE
        self._select_tab("matrix", "test1")

        ## Go back to CSD & check we are in LINEAR
        self._select_tab("csd")
        status = self._get_status()
        assert(status["selected_dataset"] == "Linear")

        ## Go back to DAG & check we are in AND
        self._select_tab("dag")
        status = self._get_status()
        assert(status["selected_dataset"] == "AND")

        ## Go back to MATRIX & check we are in EXAMPLE
        self._select_tab("matrix")
        status = self._get_status()
        assert(status["selected_dataset"] == "test1")

    def test_load_csv_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "good_csd.csv"))
        sleep(2)
        status = self._get_status()
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "good")
        assert(status["selected_input2build"] == "csd")
        assert(status["gene_names"] == ["A1", "B2", "C3", "D4"])

    def test_load_corrupt_csv_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "corrupt_csd.csv"))
        error_message = self._get_error_message()
        assert(error_message, "Your csv data can not be loaded. Make sure it only contains 0 and 1.")

    def test_load_CSD_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "csd_good.RDS"))
        sleep(2)
        status = self._get_status()
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "CSD_custom")
        assert(status["selected_input2build"] == "csd")
        assert(status["gene_names"] == ["A1", "B2", "C3", "D4"])

    def test_load_corrupt_CSD_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "csd_bad.RDS"))
        error_message = self._get_error_message()
        assert(error_message, "Your csv data can not be loaded. Make sure it only contains 0 and 1.")
    

    def test_load_DAG_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "dag_good.RDS"))
        sleep(2)
        status = self._get_status("dag")
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "DAG_custom")
        assert(status["selected_input2build"] == "dag")
        assert(status["gene_names"] == ["A1", "B2", "C3", "D4"])
    
    def test_load_corrupt_DAG_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "dag_bad.RDS"))
        sleep(2)
        error_message = self._get_error_message()
        assert(error_message, "There was a problem when checking your .rds file. Make sure it containis $type (either 'csd', 'dag', or 'matrix'), $data only with 0 and 1")
    

    def test_load_MATRIX_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "matrix_good.RDS"))
        sleep(2)
        status = self._get_status("matrix")
        assert(status["number_of_genes"] == 4)
        assert(status["selected_dataset"] == "MHN_custom")
        assert(status["selected_input2build"] == "matrix")
        assert(status["gene_names"] == ["A1", "B2", "C3", "D4"])
    
    def test_load_corrupt_MATRIX_dataset(self):
        upload = self.driver.find_element_by_css_selector(".upload_file input[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "matrix_bad.RDS"))
        sleep(2)
        error_message = self._get_error_message()
        assert(error_message, "There was a problem when checking your .rds file. Make sure it containis $type (either 'csd', 'dag', or 'matrix'), $data only with 0 and 1")
    
## TESTING CSD 
class test_csd_input(evamtools_basics):
    def test_modify_genotype_with_buttons(self):
        ## Initial status
        sleep(0.2)
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
                    "#genotype input[value=A]").click()
        sleep(1.5)
        current_freq = self.driver.find_element_by_css_selector("#genotype_freq").get_attribute("value")
        sleep(1)
        assert(int(current_freq) == 100) ## This sometimes fails with invalid literal for int() with base 10: ''

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
        sleep(1.5)
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
        sleep(0.6)

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

        ## Set freqs (remove some genotypes)
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.LEFT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys(200).perform()
        actions.send_keys(Keys.TAB).perform()
        for i in [50, 150, 0]:
            actions.send_keys(i)
            actions.send_keys(Keys.TAB).perform()

        ## Save changes
        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()
        sleep(1)

        ## Check Status
        out_table_info = self._get_table_info()
        out_genotypes = self._process_genotype_table(out_table_info)

        assert(out_genotypes["WT"] == 200)
        assert(out_genotypes["A"] == 50)
        assert(out_genotypes["AB"] == 150)
        assert(("ABC" in out_genotypes.keys()) == False)
        assert(out_genotypes["AC"] == prev_genotypes["AC"])
        assert(out_genotypes["ABCD"] == prev_genotypes["ABCD"])

    def test_change_gene_number(self):
        ## Select Linear & save data
        self.driver.set_window_size(1366, 768)
        sleep(0.2)
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

    def test_change_gene_names(self): 
        ## Select Linear & save data
        sleep(0.2)
        self._select_tab("csd", "Linear")
        sleep(0.2)
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
        self._change_gene_names(new_names)

        ## Get gene names
        status = self._get_status()
        assert(status["gene_names"], new_names)

        ## Get genotypes freqs
        genotype_info = self._get_table_info()
        genotype_info_dict = self._process_genotype_table(genotype_info)

        # Failed for no reason
        assert(genotype_info_dict == expected_genotypes_dict)

    def test_save_data_set(self):
        ## Selecting LINEAR
        sleep(0.2)
        self._select_tab("csd", "Linear")
        sleep(0.5)
        initial_genotypes = self._get_table_info()
        initial_genotypes_dict = self._process_genotype_table(initial_genotypes)

        ## Introducing modifications
        changes = {"A, B": 1000, "A": 500}

        for (genotype, freq) in changes.items():
            self._modify_genotype(genotype, freq)
            sleep(2.5)

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
        sleep(0.2)
        self._select_tab("csd", "Linear")
        final_genotypes = self._get_table_info()

        assert(final_genotypes == initial_genotypes)

# TESTING DAG
class test_dag_input(evamtools_basics):
    def test_modify_dag(self):
        ## Increasing number of genes
        self.driver.set_window_size(1366, 768)
        sleep(0.2)
        self._select_tab("dag")
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()

        # Adding edges & checking changes
        self._add_edge("Root", "B")
        self._add_edge("B", "C")
        dag_info = self._get_dag_info()
        current_dag = [
            ["Root", "A", "Single", "1"],
            ["Root", "B", "Single", "1"],
            ["B", "C", "Single", "1"],
        ]
        assert(current_dag == dag_info)

        ## Same parent and children
        self._add_edge("B", "B")
        error_message = self._get_error_message()
        sleep(0.5)
        self.driver.find_element_by_css_selector(".modal-open .modal-footer button").click()
        assert(error_message == "Both From and To options must be different")

        ## One children of root with multiple parents
        self._add_edge("D", "B")
        error_message = self._get_error_message()
        self.driver.find_element_by_css_selector(".modal-open .modal-footer button").click()
        assert(error_message == "A direct children of Root cannot have multiple parents")

        # Adding repeated relathionship
        self._add_edge("B", "C")
        error_message = self._get_error_message()
        self.driver.find_element_by_css_selector(".modal-open .modal-footer button").click()
        assert(error_message == "That edge is already present")

        # Adding bidirectional relathionship
        self._add_edge("C", "D")
        self._add_edge("D", "C")
        error_message = self._get_error_message()
        self.driver.find_element_by_css_selector(".modal-open .modal-footer button").click()
        assert(error_message == "Relationships cannot be bidirectional")

        # Adding cycle
        self._add_edge("D", "E")
        self._add_edge("E", "C")
        error_message = self._get_error_message()
        self.driver.find_element_by_css_selector(".modal-open .modal-footer button").click()
        assert(error_message == "This relationship breaks the DAG. Revise it.")

    def test_modify_dag_2(self):
        ## Increasing number of genes
        self.driver.set_window_size(1366, 768)
        sleep(0.2)
        self._select_tab("dag")
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()
        
        # Adding relationship from unlinked node
        self._add_edge("D", "E")
        current_dag = [
            ["Root", "A", "Single", "1"],
            ["Root", "D", "Single", "1"],
            ["D", "E", "Single", "1"],
        ]
        dag_info = self._get_dag_info()
        assert(current_dag == dag_info)

        ## Analysis button is enabled after sampling
        self._scroll2top()
        analysis_button = self.driver.find_element_by_css_selector("#analysis")
        assert(analysis_button.is_enabled() == False)
        
        self.driver.find_element_by_css_selector("#resample_dag").click()
        sleep(3)
        assert(analysis_button.is_enabled())

    def test_modify_gene_names(self):
        ## Increasing number of genes
        self.driver.set_window_size(1366, 768)
        sleep(0.2)
        self._select_tab("dag")
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        sleep(0.2)
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()

        # Adding edges & checking changes
        self._add_edge("Root", "B")
        self._add_edge("B", "C")
        self._add_edge("D", "E")
        
        # Changing gene names
        new_gene_names = ["A", "B1", "C", "D3"]
        expected_dag = [
            ["Root", "A", "Single", "1"],
            ["Root", "B1", "Single", "1"],
            ["B1", "C", "Single", "1"],
            ["Root", "D3", "Single", "1"],
            ["D3", "E", "Single", "1"],
        ]
        self._change_gene_names(new_gene_names)

        assert(expected_dag == self._get_dag_info())

    def test_modify_dag_from_table(self):
        ## Selecting dag AND dataset
        sleep(1)
        self._select_tab("dag", "AND")
        sleep(1)
        initial_dag = self._get_dag_info()

        ## Modify dataset
        actions = ActionChains(self.driver)
        first_cell = self.driver.find_elements_by_css_selector("#dag_table tbody td")[2]

        actions.double_click(first_cell).perform()
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.RIGHT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys("OR").perform()
        actions.send_keys(Keys.TAB).perform()
        
        actions.send_keys(4).perform()

        for _ in range(0,3):
            actions.send_keys(Keys.TAB).perform()
        actions.send_keys("bad_relationship").perform()
        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(-7).perform()

        for _ in range(0,7):
            actions.send_keys(Keys.TAB).perform()
        
        actions.send_keys("XOR").perform()
        actions.send_keys(Keys.TAB).perform()
        
        actions.send_keys(1).perform()

        for _ in range(0,3):
            actions.send_keys(Keys.TAB).perform()
        actions.send_keys("OR").perform()
        ## Saving
        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()
        sleep(2)
        ## Checking
        dag_info = self._get_dag_info()
        expected_dag_info = [
            ["Root", "A", "Single", "4"],
            ["A", "B", "Single", "2"],
            ["A", "C", "Single", "3"],
            ["B", "D", "OR", "1"],
            ["C", "D", "OR", "1"],
        ]
        assert(dag_info == expected_dag_info)

    def test_remove_node(self):
    ## Selecting dag AND dataset
        sleep(0.2)
        self._select_tab("dag", "AND")
        sleep(1)
        initial_dag = self._get_dag_info()
    
    ## Modify dataset
        actions = ActionChains(self.driver)
        first_cell = self.driver.find_elements_by_css_selector("#dag_table tbody td")[11]

        actions.double_click(first_cell).perform()
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.LEFT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys(0).perform()
        actions.send_keys(Keys.TAB).perform()
        
        for _ in range(0,3):
            actions.send_keys(Keys.TAB).perform()
        actions.send_keys(0).perform()

        for _ in range(0,4):
            actions.send_keys(Keys.TAB).perform()
        actions.send_keys(0).perform()

        ## Saving
        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()
        sleep(0.2)
        ## Checking
        dag_info = self._get_dag_info()
        expected_dag_info = [
            ["Root", "A", "Single", "1"],
            ["A", "B", "Single", "2"],
        ]
        assert(dag_info == expected_dag_info)

    def test_remove_node_2(self):
    ## Selecting dag AND dataset
        sleep(0.2)
        self._select_tab("dag", "AND")
        sleep(1)
        initial_dag = self._get_dag_info()
    
    ## Modify dataset
        actions = ActionChains(self.driver)
        first_cell = self.driver.find_elements_by_css_selector("#dag_table tbody td")[11]

        actions.double_click(first_cell)
        actions.perform()
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.LEFT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys(0).perform()
        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()
        sleep(1.5)
        ## Checking
        dag_info = self._get_dag_info()
        expected_dag_info = [
            ["Root", "A", "Single", "1"],
            ["A", "B", "Single", "2"],
            ["Root", "C", "Single", "3"],
            ["B", "D", "AND", "4"],
            ["C", "D", "AND", "4"],
        ]
        assert(dag_info == expected_dag_info)

        first_cell = self.driver.find_elements_by_css_selector("#dag_table tbody td")[15]
        actions = ActionChains(self.driver)
        actions.double_click(first_cell)
        actions.perform()
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.LEFT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys(0).perform()

        for _ in range(0,4):
            actions.send_keys(Keys.TAB).perform()
        actions.send_keys(0).perform()

        ## Saving
        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()
        sleep(0.5)
        ## Checking
        dag_info = self._get_dag_info()
        expected_dag_info = [
            ["Root", "A", "Single", "1"],
            ["A", "B", "Single", "2"],
            ["Root", "C", "Single", "3"],
        ]
        assert(dag_info == expected_dag_info)

    def test_save_data_set(self):
        ## Increasing number of genes
        self.driver.set_window_size(1366, 768)
        sleep(0.2)
        self._select_tab("dag")
        old_dag = self._get_dag_info()
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()
        sleep(0.5)

        # Adding edges & checking changes
        self._add_edge("Root", "B")
        self._add_edge("B", "C")
        self._add_edge("D", "E")

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
        assert(status["selected_input2build"] == "dag")

        new_dag = self._get_dag_info()
        expected_dag = [
            ["Root", "A", "Single", "1"],
            ["Root", "B", "Single", "1"],
            ["B", "C", "Single", "1"],
            ["Root", "D", "Single", "1"],
            ["D", "E", "Single", "1"],
        ]
        assert(expected_dag == new_dag)

        ## Check dataset restoration
        self._select_tab("dag", "User")
        user_dag = self._get_dag_info()

        assert(old_dag == user_dag)

## TESTING MATRIX
class test_matrix_input(evamtools_basics):
    def test_change_theta(self):
        ## Selecting dag AND dataset
        sleep(1)
        self._select_tab("matrix", "test1")
        sleep(1)

        analysis_button = self.driver.find_element_by_css_selector("#analysis")
        assert(analysis_button.is_enabled() == False)
        ## Modify dataset
        actions = ActionChains(self.driver)
        first_cell = self.driver.find_elements_by_css_selector("#thetas_table tbody td")[2]


        actions.double_click(first_cell).perform()
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.LEFT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys("2").perform()
        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(4).perform()

        for _ in range(0,3):
            actions.send_keys(Keys.TAB).perform()
        actions.send_keys("0").perform()
        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(-7).perform()

        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()
        sleep(1)

        expected_theta_table = [
            ["A1", "B2", "C3"],
            ["A1", "-0.4", "2", "4"],
            ["B2", "-1.04", "0", "-7"],
            ["C3", "-2.53", "-0.76", "0.23"],
        ]

        theta_table = self._get_theta_table()
        assert(theta_table == expected_theta_table)

        assert(analysis_button.is_enabled())

    def test_change_names(self):
        self._select_tab("matrix")
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(50, 0).release().perform()
        sleep(0.5)

        new_gene_names = ["A1", "B", "C3"]
        self._change_gene_names(new_gene_names)
        theta_table = self._get_theta_table()

        expected_theta_table = [
            ["A1", "B", "C3", "D"],
            ["A1", "0", "0", "0", "0"],
            ["B", "0", "0", "0", "0"],
            ["C3", "0", "0", "0", "0"],
            ["D", "0", "0", "0", "0"],
        ]
        assert(theta_table == expected_theta_table)

    def test_change_gene_number(self):
        self._select_tab("matrix")
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()
        sleep(0.5)

        status = self._get_status("matrix")

        assert(status["gene_names"] == ["A", "B", "C", "D", "E", "F"])

    def test_save_data_set(self):
        self._select_tab("matrix", "test1")
        old_thetas = self._get_theta_table()

        # Changing gene number
        if self.driver.get_window_size()["width"] > 1500:
            self.driver.set_window_size(1366, 768)
        sleep(0.2)
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(60, 0).release().perform()
        sleep(0.5)
  
        ## Modify dataset
        actions = ActionChains(self.driver)
        first_cell = self.driver.find_elements_by_css_selector("#thetas_table tbody td")[2]

        actions.double_click(first_cell).perform()
        sleep(0.4)
        actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.LEFT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).send_keys(2).perform()
        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(4).perform()

        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(3).perform()

        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(2).perform()

        for _ in range(0,3):
            actions.send_keys(Keys.TAB).perform()
        actions.send_keys("0").perform()
        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(-7).perform()

        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(1).perform()

        actions.send_keys(Keys.TAB).perform()
        actions.send_keys(-1).perform()

        actions.key_down(Keys.CONTROL).send_keys(Keys.ENTER).key_up(Keys.CONTROL).perform()
        sleep(1)

        # Changing genes names
        new_gene_names = ["A1", "B", "C3", "D", "E5"]
        self._change_gene_names(new_gene_names)

        expected_theta_table = [
            ["A1", "B", "C3", "D", "E5"],
            ["A1", "-0.4", "2", "4", "3", "2"],
            ["B", "-1.04", "0", "-7", "1", "-1"],
            ["C3", "-2.53", "-0.76", "0.23", "0", "0"],
            ["D", "0", "0", "0", "0", "0"],
            ["E5", "0", "0", "0", "0", "0"],
        ]

        theta_table = self._get_theta_table()
        assert(theta_table == expected_theta_table)

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
        status = self._get_status("matrix")
        assert(status["selected_dataset"] == new_dataset_name)
        assert(status["selected_input2build"] == "matrix")

        theta_table = self._get_theta_table()
        assert(theta_table == expected_theta_table)

        ## Check dataset restoration
        self._select_tab("matrix", "test1")
        sleep(1)
        theta_table = self._get_theta_table()
        sleep(0.1)
        assert(theta_table == old_thetas) ## Failed then OK

        self._select_tab("dag")
        sleep(1)
        self._select_tab("matrix", new_dataset_name)
        sleep(1)
        theta_table = self._get_theta_table()
        assert(theta_table == expected_theta_table)

## TESTING RESULTS
class test_results(evamtools_basics):
    def _check_tabular_data(self, input2build, do_mccbn = False):
        available_cpms = self.driver.find_elements_by_css_selector("#cpm2show .checkbox input")
        available_cpms = [i.get_attribute("value") for i in available_cpms]
        available_cpms.remove("Source")
        # if not do_mccbn: 
        #     available_cpms.remove("MCCBN")

        not_ot = available_cpms.copy()
        not_ot.remove("OT")
        # not_ot.remove("DBN")

        expected_tabular = not_ot.copy()

        if input2build != "csd":
            expected_tabular.append("Source")
        expected_tabular_from_sims = ["From", "To", *expected_tabular]

        expected_lambdas = available_cpms.copy()
        if input2build == "dag":
            expected_lambdas.append("Source")
        expected_lambdas.remove("MHN")
        expected_lambdas = ["Gene", *expected_lambdas]

        expected_td_trans_mat = expected_tabular_from_sims.copy()
        if input2build != "csd":
            expected_td_trans_mat.remove("Source")
            ## FIXME : check next

        tabular_types = {
            # "f_graph": expected_tabular_from_sims,
            "trans_rate_mat": expected_tabular_from_sims,            
            "genotype_transitions": expected_tabular_from_sims,
            "freqs": ["Genotype", "OT", *expected_tabular],
            "trans_mat": ["From", "To", "OT", *expected_tabular],
            "lambdas": expected_lambdas,
            "td_trans_mat": expected_td_trans_mat,
        }

        for k in tabular_types.keys():
            self.driver.find_element_by_css_selector(f"#data2table input[value='{k}']").click()
            sleep(0.4)
            table_info = self._get_cpm_table()
            print(k)
            print("what I want", tabular_types[k])
            print("what I see ", table_info[0])
            assert(set(table_info[0]) == set(tabular_types[k]))

    def _basic_results_test(self, data):
        ## Check status
        active_cpms = self.driver.find_elements_by_css_selector("#cpm2show .checkbox input[checked=checked]")
        active_cpms = [i.get_attribute("value") for i in active_cpms]
        for i in active_cpms:
            plot_sims_1 = self.driver.find_elements_by_css_selector(f"#plot_sims_{i}")
            assert(len(plot_sims_1) == 1)
            plot_sims_2 = self.driver.find_elements_by_css_selector(f"#plot_sims_{i}")
            assert(len(plot_sims_2) == 1)
        
        ## Check name of downloaded file
        download_button = self.driver.find_element_by_css_selector("#download_cpm").click()
        sleep(1)
        base_folder = os.getcwd()
        os.chdir(os.path.join(os.path.expanduser("~"), "Downloads"))
        last_created_file = sorted(filter(os.path.isfile, os.listdir('.')), key=os.path.getmtime)[-1]
        assert(last_created_file[-4:] == ".RDS")
        os.remove(last_created_file)
        
        ## Check wether "source" should be dissabled
        source_option = self.driver.find_element_by_css_selector("#cpm2show input[value=Source]")
        assert(source_option.is_enabled() == (data["tab"] != "csd" ))

        ## Check wether "MCCBN" should be dissabled
        # source_option = self.driver.find_element_by_css_selector("#cpm2show input[value=MCCBN]")
        # assert(source_option.is_enabled() == data["mccbn"])

        ## Check functionality of 'Modify data' button
        self.driver.find_element_by_css_selector("#modify_data").click()
        sleep(1.5)

        status = self._get_status(data["tab"])
        assert(status["number_of_genes"] == data["number_of_genes"])
        assert(status["selected_dataset"] == data["name"])
        assert(status["selected_input2build"] == data["tab"])
        assert(status["gene_names"] == data["gene_names"])
        os.chdir(base_folder)

    def _get_cpm_table(self):
        row_table = self.driver.find_elements_by_css_selector("#cpm_freqs tr")
        row_text = [ i.text.replace(", ", ",").split(" ") for i in row_table ]
        return row_text

    def test_results_on_load(self):
        ## Switch to results tab
        self._go_to("result_viewer")
        
        ## Placeholder text to ask for running something
        # sleep(1)
        # placeholder_text = self.driver.find_element_by_css_selector("#sims2 h3").text
        # assert(placeholder_text == "There are not results to show yet. Go to the input tab, select a dataset and hit the 'Run evamtools!' button")

        # download_button = self.driver.find_elements_by_css_selector("#download_cpm[disabled=disabled]")
        # assert(len(download_button) == 1)
        # assert(download_button[0].text == 'Download!')

        ## Uploading data
        upload = self.driver.find_element_by_css_selector("input#output_cpms[type=file]")
        upload.send_keys(os.path.join(os.getcwd(), "AND_test_cpm.RDS")) 
        sleep(1)

        ## Status changes after uploading some data

        download_button = self.driver.find_elements_by_css_selector("#download_cpm")
        assert(len(download_button) == 1)
        assert(download_button[0].text == 'Download!')
        download_button = self.driver.find_elements_by_css_selector("#download_cpm[disabled=disabled]")
        assert(len(download_button) == 0)

        expected_data = {"tab": "csd", "gene_names": ["A", "B", "C", "D"],
            "name": "AND_new", "number_of_genes": 4, "mccbn": False}
        self._basic_results_test(expected_data)     

    def test_results_on_CSD(self):
        ## Load data
        self._select_tab("csd", "AND")
        sleep(0.5)

        ## Modify data
        self._modify_genotype("A", 500)
        self._modify_genotype("A, B", 50)

        ## Save data 
        new_dataset_name = "AND_new"
        dataset_name = self.driver.find_element_by_css_selector("input#dataset_name")
        dataset_name.clear()
        dataset_name.send_keys(new_dataset_name)
        sleep(0.5)
        save_button = self.driver.find_element_by_id("save_csd_data")
        save_button.click()
        sleep(1)

        ## Run analysis
        self._scroll2top()
        self.driver.find_element_by_css_selector("#analysis").click()
        sleep(1)
        self._select_tab("csd", "User")
        sleep(10)

        ## Check new status
        self._go_to("result_viewer")
        self._select_result(new_dataset_name)
        
        expected_data = {"tab": "csd", "gene_names": ["A", "B", "C", "D"],
            "name": new_dataset_name, "number_of_genes": 4}

        self._basic_results_test(expected_data) 
        self._go_to("result_viewer")
        self._check_tabular_data(expected_data["tab"])

    def test_results_on_DAG(self):
        ## Load data
        self._select_tab("dag")
        sleep(0.5)

        ## Modify data
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()
        sleep(0.5)

        self._add_edge("Root", "B")
        self._add_edge("B", "C")
        self._add_edge("C", "D")
        new_gene_names = ["A", "B1", "C", "D3"]
        self._change_gene_names(new_gene_names)

        self.driver.find_element_by_css_selector("#resample_dag").click()
        sleep(5)

        ## Save data 
        new_dataset_name = "DAG_new"
        dataset_name = self.driver.find_element_by_css_selector("input#dataset_name")
        dataset_name.clear()
        dataset_name.send_keys(new_dataset_name)
        sleep(0.5)
        save_button = self.driver.find_element_by_id("save_csd_data")
        save_button.click()
        sleep(1)
        
        ## Run analysis
        self._scroll2top()
        self.driver.find_element_by_css_selector("#analysis").click()
        sleep(1)
        self._select_tab("csd", "User")
        sleep(10)

        ## Check new status
        self._go_to("result_viewer")
        self._select_result(new_dataset_name)
        
        expected_data = {"tab": "dag", "gene_names": ["A", "B1", "C", "D3"],
            "name": new_dataset_name, "number_of_genes": 4, "mccbn": False}
        self._basic_results_test(expected_data) 
        self._go_to("result_viewer")
        self._check_tabular_data(expected_data["tab"], expected_data["mccbn"])

    def test_results_on_matrix(self):
        ## Load data
        self._select_tab("matrix", "test1")
        sleep(0.5)

        ## Sampling
        self.driver.find_element_by_css_selector("#resample_mhn").click()
        sleep(1)

        ## Save data 
        new_dataset_name = "MHN_new"
        dataset_name = self.driver.find_element_by_css_selector("input#dataset_name")
        dataset_name.clear()
        dataset_name.send_keys(new_dataset_name)
        sleep(0.5)
        save_button = self.driver.find_element_by_id("save_csd_data")
        save_button.click()
        sleep(1)

        ## Run analysis
        self._scroll2top()
        self.driver.find_element_by_css_selector("#analysis").click()
        sleep(1)
        self._select_tab("csd", "User")
        sleep(10)

        ## Check new status
        self._go_to("result_viewer")
        self._select_result(new_dataset_name)
        
        expected_data = {"tab": "matrix", "gene_names": ["A", "B", "C"],
            "name": new_dataset_name, "number_of_genes": 3, "mccbn": False}
        self._basic_results_test(expected_data) ## Sometimes fails for no apparent reason
        self._go_to("result_viewer")
        self._check_tabular_data(expected_data["tab"], expected_data["mccbn"])
    
    # def test_running_MCCBN(self):
    #     pdb.set_trace()
    #     self.driver.find_element_by_css_selector("#advanced_options").click()
    #     sleep(1)
    #     self.driver.find_element_by_css_selector("#more_cpms input[value='mccbn']").click()
        # ## Load data
        # self._select_tab("matrix", "test1")
        # sleep(0.5)

        # ## Sampling
        # self.driver.find_element_by_css_selector("#resample_mhn").click()
        # sleep(1)

        # ## Save data 
        # new_dataset_name = "MHN_new"
        # dataset_name = self.driver.find_element_by_css_selector("input#dataset_name")
        # dataset_name.clear()
        # dataset_name.send_keys(new_dataset_name)
        # sleep(0.5)
        # save_button = self.driver.find_element_by_id("save_csd_data")
        # save_button.click()
        # sleep(1)

        # self.driver.find_element_by_tag_name('body').send_keys(Keys.CONTROL + Keys.HOME)
        # ## Run MCCBN
        # self.driver.find_element_by_css_selector("#advanced_options").click()
        # sleep(1)
        # self.driver.find_element_by_css_selector("#more_cpms input[value='mccbn']").click()
        # sleep(0.5)
        # self.driver.find_element_by_css_selector(".modal-open .modal-footer button").click()

        # ## Run analysis
        # self.driver.find_element_by_css_selector("#analysis").click()
        # sleep(1)
        # self._select_tab("csd", "User")
        # sleep(10)

        # ## Check new status
        # self._go_to("result_viewer")
        # self._select_result(new_dataset_name)
        
        # expected_data = {"tab": "matrix", "gene_names": ["A", "B", "C"],
        #     "name": new_dataset_name, "number_of_genes": 3, "mccbn": True}

        # self._basic_results_test(expected_data) 
        # self._go_to("result_viewer")
        # self._check_tabular_data(expected_data["tab"], expected_data["mccbn"])

if __name__ == "__main__":
    unittest.main(failfast = True)
