from imagepipe.raw_functions import split_and_trim
current_location= "D:\\AKN_TIFF 8BIT"
main_root= "D:\\AKN_TIFF 8BIT\\hello"
prefix = split_and_trim(current_location, main_root)
print main_root[0:]