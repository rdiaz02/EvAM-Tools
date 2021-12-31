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
        return self.driver.find_element_by_css_selector("#shiny-modal .modal-body div").text
    
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

    def _select_tab(self, input_type, name_dataset = None):
        csd_tab = self.driver.find_element_by_css_selector(f"#input2build .radio input[value={input_type}]")
        csd_tab.click()
        sleep(0.5)
        if(name_dataset):
            dataset_tab = self.driver.find_element_by_css_selector(f"#select_csd .radio input[value={name_dataset}]")
            dataset_tab.click()
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
        error_message = self._get_error_message()
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
    
    ## TESTING CSD 
    def test_modiying_genotype_with_buttons(self):
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
                    "#genotype input[value=A]")
        input_gene.find_element_by_xpath('..').click()
        current_freq = self.driver.find_element_by_css_selector("#genotype_freq").get_attribute("value")
        print(current_freq)
        sleep(1)
        if(current_freq != "100"):
            pdb.set_trace()
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
        sleep(0.2)
        self._select_tab("csd", "AND")
        sleep(0.2)

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
        #actions.key_down(Keys.CONTROL).key_down(Keys.SHIFT).send_keys(Keys.RIGHT).key_up(Keys.CONTROL).key_up(Keys.SHIFT).perform()

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

    def test_change_gene_number_CSD(self):
        ## Select Linear & save data
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

    def test_change_gene_names_CSD(self):
        ## Select Linear & save data
        sleep(0.2)
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
        self._change_gene_names(new_names)

        ## Get gene names
        status = self._get_status()
        assert(status["gene_names"], new_names)

        ## Get genotypes freqs
        genotype_info = self._get_table_info()
        genotype_info_dict = self._process_genotype_table(genotype_info)

        assert(genotype_info_dict == expected_genotypes_dict)

    def test_save_data_set_CSD(self):
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
        self._select_tab("csd", "Linear")
        final_genotypes = self._get_table_info()

        assert(final_genotypes == initial_genotypes)

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

    # TESTING DAG
    def test_modifying_dag(self):
        ## Increasing number of genes
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

    def test_more_dag_modifations(self):
        ## Increasing number of genes
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
        analysis_button = self.driver.find_element_by_css_selector("#analysis")
        assert(analysis_button.is_enabled() == False)
        
        self.driver.find_element_by_css_selector("#resample_dag").click()
        sleep(3)
        assert(analysis_button.is_enabled())

    def test_change_names_DAG(self):
        ## Increasing number of genes
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

    def test_modify_dag_value_from_table(self):
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
        # if(dag_info != expected_dag_info):
            # pdb.set_trace()
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

    def test_more_remove_node(self):
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

    def test_save_data_set_DAG(self):
        ## Increasing number of genes
        self._select_tab("dag")
        old_dag = self._get_dag_info()
        slider_input = self.driver.find_element_by_css_selector("#genes_number span.irs-handle")
        move = ActionChains(self.driver)
        move.click_and_hold(slider_input).move_by_offset(100, 0).release().perform()

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
