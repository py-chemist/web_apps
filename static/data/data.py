def get_checkboxes_data():

    title = ["Wrap the code into \chemfig{...}", "Assign number to each atom except hydrogen", "Show nicer double and triple bonds",
             "Show circle in aromatic compounds instead of double bonds", "Show carbon atoms as elements", "Show methyl group as elements",
             "Flip the structure horizontally", "Flop the strcture vertically"]
    value = ["-w", "-n", "-f", "-o", "-c", "-m", "-p", "-q"]
    text = ["chemfig", "atom-numbers", "fancy bonds", "aromatic", "show carbon", "show methyl", "flip", "flop"]
    return zip(title, value, text)
    

def get_menu_links():

    return ["Home", "Links", "About"]

def hydrogens_options():

    return [("keep", "keep"), ("add", "add"), ("delete", "delete")]    
    
        
