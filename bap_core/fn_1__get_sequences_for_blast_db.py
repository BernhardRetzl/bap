def function_4_1(root_dir, family, get_only_values, database, database_name, my_min_values, my_max_values):

    def get_sequences_for_blast_db():
        to_write = list()
        file_list = glob.glob(database)
        for file in file_list:
            print(file)
            handler = open(file)
            info = handler.readline().strip()
            sequence = handler.readline().strip()
            while info != '':
                spacing = info.split(' ')[-1].split('-')[1:-1]
                spacing = np.array([int(i) for i in spacing])
                spacing_excess_length = len(spacing) - 4
                for i in range(spacing_excess_length):
                    new_spacing = spacing[0 + i:5 + i]
                    if all(new_spacing >= new_minimum_values) and all(new_spacing <= new_maximum_values):
                        to_write.append(info)
                        to_write.append(sequence)
                        break
                info = handler.readline().strip()
                sequence = handler.readline().strip()

        with open(root_dir+'Data/'+family+'/_4_1/'+database_name, 'wt') as out_file:
            for thing in to_write:
                out_file.write(thing + '\n')

    def get_recommended_sequences():
        with open(root_dir+'Data/'+family+'/_4_0/queries.txt') as in_file:
            final_dict = dict()
            for line in in_file:
                if line.startswith('>'):
                    continue
                elif line.count('C') == 6:
                    for number, part in enumerate(line.split('C')):
                        if number in final_dict:
                            final_dict[number].append(len(part))
                        else:
                            final_dict[number] = [len(part),]
            minimum_list = list()
            maximum_list = list()
            print(final_dict)
            for item in final_dict:
                if item == 0 or item == 6:
                    continue
                else:
                    print('Loop : '+str(item)+'\n'+'The minimum value is: '+str(min(final_dict[item]))+
                          '\n'+'The maximum value is: '+str(max(final_dict[item])))
                    minimum_list.append(min(final_dict[item]))
                    maximum_list.append(max(final_dict[item]))

    get_recommended_sequences()
    if not get_only_values:
        new_minimum_values = np.array(my_min_values)
        new_maximum_values = np.array(my_max_values)
        get_sequences_for_blast_db()