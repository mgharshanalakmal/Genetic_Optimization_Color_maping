import random
import numpy as np

class gen_optimization:
    
    def __init__(self, num_rows, num_cols, num_colors, chromosome_count, best_parent) -> None:
        """Initializing the class.

        Args:
            num_rows (int): Number of rows in the grid. ex:5 
            num_cols (int): Number of columns in the grid. ex:5
            num_colors (int): Number of colors to be used when coloring cells. ex:7
            chromosome_count (int): Number of chromosomes in the population ex:1000
            best_parent (int): Number of best parents to be selected from the selection process. ex:100 
        """        
        
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.num_colors = num_colors
        self.chromosome_count = chromosome_count
        self.best_parent = best_parent
        self.adjacent_values = self.get_adjacent_cells(self.num_rows, self.num_cols)
        self.chromosome_list = self.create_random_solution(self.num_rows, self.num_cols, self.num_colors)
    
    
    def get_adjacent_cells(self, num_rows, num_cols):
        """Find all adjacent cell for the each cell of the grid and stored values in a dictionary to be used while the selection process 
        to find the fitness of each chromosome.  

        Args:
            num_rows (int): Number of rows in the grid
            num_cols (int): Number of columns in the grid

        Returns:
            Dictionary: Dictionary that contains list of adjacent cells of each cell. Cell numbers are the keys in the dictionary. Values stored as lists.
        """        
    
        num_cells = num_rows * num_cols
        values = list(range(num_cells))
        grid = [values[i:i+num_cols] for i in range(0, num_cells, num_cols)]
        
        
        adjacent_cells = {}
        
        for row in range(num_rows):
            for col in range(num_cols):
                
                cell_num = row * num_cols + col
                neighbors = []
                    
                #First Row
                if row == 0:
                    #Lower Cell
                    neighbors.append((row + 1) * num_cols + col)
                    
                    #Upper Left Corner
                    if col == 0:
                            #Right cell
                            neighbors.append(row * num_cols + col + 1)
                            #Lower Right
                            neighbors.append((row + 1) * num_cols + col + 1)
                        
                    #Upper Right Corner
                    elif col == (num_cols - 1):
                            #Left cell
                            neighbors.append(row * num_cols + col - 1)
                            #Lower Left
                            neighbors.append((row + 1) * num_cols + col - 1)
                            
                    else:
                            
                            #Right cell
                            neighbors.append(row * num_cols + col + 1)
                            #Left cell
                            neighbors.append(row * num_cols + col - 1)
                            #Lower Right
                            neighbors.append((row + 1) * num_cols + col + 1)
                            #lower Left
                            neighbors.append((row + 1) * num_cols + col - 1)
                
                #Last Row           
                elif row == num_rows - 1:
                    #Upper Cell
                    neighbors.append((row - 1) * num_cols + col)
                    
                    #Lower Left Corner
                    if col == 0:
                            #Right cell
                            neighbors.append(row * num_cols + col + 1)
                            #Upper Right
                            neighbors.append((row - 1) * num_cols + col + 1)
                            
                    #Lower Right Corner
                    elif col == (num_cols - 1):
                            #Left cell
                            neighbors.append(row * num_cols + col - 1)
                            #Upper Left
                            neighbors.append((row - 1) * num_cols + col - 1)
                            
                    else:
                            #Right cell
                            neighbors.append(row * num_cols + col + 1)
                            #Upper Right
                            neighbors.append((row - 1) * num_cols + col + 1)
                            #Left cell
                            neighbors.append(row * num_cols + col - 1)
                            #Upper Left
                            neighbors.append((row - 1) * num_cols + col - 1)
                            
                #Left Column
                elif col == 0 and row > 0 and row < (num_cols - 1):
                    #Upper cell
                    neighbors.append((row - 1) * num_cols + col)
                    #Upper Right
                    neighbors.append((row - 1) * num_cols + col + 1)
                    #Right cell
                    neighbors.append(row * num_cols + col + 1)
                    #Lower Right
                    neighbors.append((row + 1) * num_cols + col + 1)
                    #Lower
                    neighbors.append((row + 1) * num_cols + col)
                    
                #Right Columnn
                elif col == (num_cols - 1) and row > 0 and row < (num_cols - 1):
                    #Upper cell
                    neighbors.append((row - 1) * num_cols + col)
                    #Upper Left
                    neighbors.append((row - 1) * num_cols + col - 1)
                    #Left cell
                    neighbors.append(row * num_cols + col - 1)
                    #Lower Left
                    neighbors.append((row + 1) * num_cols + col - 1)
                    #Lower cell
                    neighbors.append((row + 1) * num_cols + col)
                    
                                    
                else:
                    #Upper cell
                    neighbors.append((row - 1) * num_cols + col)
                    #Upper Right
                    neighbors.append((row - 1) * num_cols + col + 1)
                    #Right cell
                    neighbors.append(row * num_cols + col + 1)
                    #Lower Right
                    neighbors.append((row + 1) * num_cols + col + 1)
                    #Lower
                    neighbors.append((row + 1) * num_cols + col)
                    #Lower Left
                    neighbors.append((row + 1) * num_cols + col - 1)
                    #Left cell
                    neighbors.append(row * num_cols + col - 1)
                    #Upper Left
                    neighbors.append((row - 1) * num_cols + col - 1)
                
                            
                    
                adjacent_cells[cell_num] = neighbors    

        return adjacent_cells
    
    
    def create_random_solution(self, num_rows, num_cols, num_colors):
        """Create random solutions for the grid.  

        Args:
            num_rows (int): Number of rows in the grid.
            num_cols (int): Number of columns in the grid.
            num_colors (int): Number of colors to use.

        Returns:
            list: A list size of chromosome_count, that contains lists with size num_rows*num_cols to represent the each cell of the grid. 
            Each cell is assigned by a color from no of colors randomly. Output is a list that contains lists. 
        """        
        
        chromosome_list = []
        for i in range(self.chromosome_count):        
            chromosome_list.append([random.randint(0, num_colors - 1) for _ in range(num_rows * num_cols)])
            
        return chromosome_list
    
    
    def repeated_count(self, adjacent_values, color_list):
        """This is the function to use in the selction process. Here each cell color is selected and compare it with colors of adjacent cells. If the colors of the adjacent cell list 
        contain the same color of the cell, it will count as repeated colors. And through a loop repeated values will be aggregated. This is the values that use in the seleaction process
        as the value of the fitness of the chromosome.

        Args:
            adjacent_values (dictionary): Dictionary that contains adjacent cells for all the cells of the grid. This is the return of the 'get_adjacent_cells' function.
            color_list (list): List size of rows*columns that contains colors for the cells of the grid.

        Returns:
            int: number of adjacent cells with same color.
        """        
        repeated_count = 0
        
        for value in adjacent_values:
            cell_color = color_list[value]
            adjacent_cell_colors = [color_list[i] for i in adjacent_values[value]]
            repeated_color_count = adjacent_cell_colors.count(cell_color)
            repeated_count += repeated_color_count
        
        return repeated_count
    
    
    def setting_fittness(self, chromosome_list):
        """Getting repeat color cell count of each chromosome in the population generated using the function 'create_random_solutions'. 

        Args:
            chromosome_list (list): A list that contains list of grid values, generated by 'create_random_solutions' function.

        Returns:
            dictionary: Dictionary that contains chromosome, repeated cells in that chromosome, and the number of unique colors in the list of cells. 
        """        
        results_dict = {}
        
        for i in range(len(chromosome_list)):
            
            color_ist = chromosome_list[i]
            repeats = self.repeated_count(self.adjacent_values, color_ist)
            no_of_colors = len(set(color_ist))
            
            results_dict[i] = {
                'Chromosome': color_ist,
                'Repeat Count': repeats,
                'No of Colors': no_of_colors
            }
        
        return results_dict
    
    
    def select_best_chromosomes(self, results_dictionary):
        """Funtion that evaluate all the chromosomes of the population. Number of parents (selected chromosomes), are determined by the value entered as 'best_parent'
        at the class intiation.

        Args:
            results_dictionary (dictionary): Results dictionary created using 'setting_fittness' function. 

        Returns:
            list: list that contains dictionary values which are with the lowest repeated count.
        """        
        
        sorted_dict = sorted(results_dictionary.items(), key=lambda x: x[1]['Repeat Count'])
        lowest_keys = [x for x in sorted_dict][:self.best_parent]
        lowest_keys_index = [x[0] for x in lowest_keys]
        
        return lowest_keys
    
    
    def mutation(self, child):
        """Slightly mutate the child generated by crossing parents. Here it randomly select one value to change in the child and change that value to a randomly generated value.

        Args:
            child (list): list of of cell values genearted by crossing two selected chromosome.

        Returns:
            list: Sligtly mutated version of the child.
        """        
        
        mutation_gene = random.randint(0, len(child) - 1)
        mutated_gene_value = random.randint(0, self.num_colors - 1)
        child[mutation_gene] = mutated_gene_value
        
        return child
    
    
    def crossover(self, lowest_key):
        """Create offsping using selected parent chromosome. Select random two chromosome from selected chromosome and combine two to create new two child lists.

        Args:
            lowest_key (dictionary): dictionary that contained selected chromosome from 'select_best_chromosomes' function.

        Returns:
            list: List that contains list of cell values for the grid.
        """        
        
        new_generation = []
        index_range = np.arange(self.best_parent).tolist()
        
        for i in range(int(self.chromosome_count/2)):
            
            comb = random.sample(index_range, 2)
            chrom_combination = [lowest_key[x][1]['Chromosome'] for x in comb]
            parent_1, parent_2 = chrom_combination
            
            cross_point = random.randint(1, len(parent_1) - 2)
            child_1 = parent_1[:cross_point] + parent_2[cross_point:]
            child_2 = parent_2[:cross_point] + parent_1[cross_point:]
            
            child_1 = self.mutation(child_1)
            child_2 = self.mutation(child_2)
            new_generation.append(child_1)
            new_generation.append(child_2)
            
        return new_generation 
    
    
    def generation_score(self, generation_best_chromosomes):
        """Scoring the total generation based on the minimum number of colors and minimum number of total repeated cells.

        Args:
            generation_best_chromosomes (dictionaty): dictionary that contains best chromosomes selected through the selected process. 

        Returns:
            int: total repeated cells in the whole generation.
            int: minimum unique colors used.
        """        
        
        repeat_list = []
        color_list = []
        
        for i in range(len(generation_best_chromosomes)):
            repeat_list.append(generation_best_chromosomes[i][1]['Repeat Count'])
            color_list.append(generation_best_chromosomes[i][1]['No of Colors'])
        
        repeat_array = np.array(repeat_list)
        color_array = np.array(color_list)
        
        return repeat_array.sum(), color_array.min()
    
    
    def optimized_list(self, colors):
        """Part where all processes of genetic algorithm connected. Randomly generated chromosome list take as the initial values. Best chromosomes set to None, generation number to zero 
        and repeat vlaue to 10000. Then the while loop is created to run the loop until the repeated value become zero. Which means until the schomosome genes do to have repeated color values.
        Then the generation number is limited to 100 to break when there is no possible solution.
        Strongest chromosomes are selected from the selected process and new generation is created by the crossing. Then the new generation replace the intial shromosome ppulation. And loop 
        over and over till repeated values for the all selected parent lists become zero. 

        Args:
            colors (int): number of colors to use

        Returns:
            _type_: _description_
        """        
        
        chromosome_list = self.chromosome_list
        best_chromosomes = None
        color_count = 0
        repeat_value = 10000
        generation_number = 0
        
        while repeat_value != 0:
            if generation_number < 100:
                fitness_results = self.setting_fittness(chromosome_list)
                best_chromosomes = self.select_best_chromosomes(fitness_results)
                new_geneation = self.crossover(best_chromosomes)
                repeat_value = self.generation_score(best_chromosomes)[0]
                color_count = self.generation_score(best_chromosomes)[1]
                chromosome_list = new_geneation
                generation_number += 1
                
            else:
                
                best_chromosomes = None
                color_count = None
                print(f'No optimal solution at {colors} colors.....')
                break
            
        return best_chromosomes, generation_number
    
    



def main(rows, columns, no_of_color, chromosome_count, best_parents):
    """Loop over number of color options to check for optimal solution. If the number of colors do not return optimal solution those will be ignored. Color values with optimal solution
    will return number of colors, generation number and an example of solved grid selected randomly from selected chromosomes.

    Args:
        rows (_type_): _description_
        columns (_type_): _description_
        no_of_color (_type_): _description_
        chromosome_count (_type_): _description_
        best_parents (_type_): _description_

    Returns:
        _type_: _description_
    """    
    
    best_chromosomes = {}
    for i in range(1, no_of_color + 1):
        dd = gen_optimization(rows, columns, i, chromosome_count, best_parents)
        best_chrom, gen_number = dd.optimized_list(i)
        
        if best_chrom is None:
            pass
        
        else:
            optimal_color_list = []
            for j in range(len(best_chrom)):
                optimal_color = best_chrom[j][1]['Chromosome']
                optimal_color_list.append(optimal_color)
                
            best_chromosomes[i] = {
                'Color Lists': optimal_color_list,
                'Generation Number': gen_number
            } 


    def grid(row_num, col_num, list):
        grid = []
        for i in range(row_num):
            row = []
            for j in range(col_num):
                element = list[i * row_num + j]
                row.append(element)
            grid.append(row)
        
        return grid


    for key in best_chromosomes:
        colors = key
        gen = best_chromosomes[key]['Generation Number']
        rand_list_index = random.randint(0, len(best_chromosomes[key]['Color Lists']) - 1)
        rand_list = best_chromosomes[key]['Color Lists'][rand_list_index]
        result_grid = grid(rows, columns, rand_list)
        
        print()
        print(f'Optimal solution for {colors} colors achieved at {gen} generations.....')
        print()
        print('Example Result Grid')
        print()
        for row in result_grid:
            print(row)
        print()
        
        
if __name__ == '__main__':
    main(10,10, 8, 1000, 100)