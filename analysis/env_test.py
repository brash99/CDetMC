print("Environment:")
print("-----------------")

sys_values = True
if (sys_values):
    import sys
    print("sys.executable: ", sys.executable)
    print("-----------------")
    from pprint import pprint as p
    p(sys.path)

print("-----------------")