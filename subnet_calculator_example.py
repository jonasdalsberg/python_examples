### Script written for a friend who works in IT
### Script written by Jonas Dalsberg Jørgensen in March 2022


# This is a subnetting calculator program
# sub_cal_v3.py - sub
import ipaddress


def print_network_info(network_addr):
    number_of_hosts = len(list(network_addr.hosts()))
    print(f"network Address: {network_addr.with_prefixlen}\n"
          f"Network mask: {network_addr.netmask}\n"
          f"Broadcast address: {network_addr.broadcast_address}\n"
          f"Number of Host: {number_of_hosts}\n"
          f"Number of Addresses: {network_addr.num_addresses}\n")


# splits the original network into number of subnets needed
def network_based_subnetting(network_addr, number_of_subnets):
    numb_of_subnets_in_network = 1

    subnetting = list(network_addr.subnets(prefixlen_diff=numb_of_subnets_in_network))

    while len(subnetting) < int(number_of_subnets):
        subnetting = list(network_addr.subnets(prefixlen_diff=numb_of_subnets_in_network))
        numb_of_subnets_in_network = numb_of_subnets_in_network + 1

    for subnet in subnetting:
        print_network_info(subnet)
    print(f"Number of subnets: {len(subnetting)}\n")


# comparing number of host in netmask with host needed (finds smallest)
def host_based_subnetting(host_needed):
    # starts with /30 that has 2 hosts
    netmask = 30
    host = 2

    while host <= host_needed:
        host = host + host + 2
        netmask = netmask - 1

    print(f"\nhost you need: {host_needed}\n"
          f"smallest netmask: /{netmask}\n"
          f"number of hosts: {host}\n")


while True:

    options = input("What do you want to do(1-4):\n"
                    "1: See information about network\n"
                    "2: Subnetting a network based on number of subnets needed\n"
                    "3: find smallest netmask based on number of hosts needed\n"
                    "4: exit program\n")

    # So I can reuse the variable "network"
    if options == "1" or options == "2":
        # Convert User input to a ipv4 network object
        # strick=False enables the user to write a host or network address
        user_network = input("\nWrite a network/ip with Prefix or netmask (addr/xx or addr/xxx.xxx.xxx.xxx):\n")
        network = ipaddress.IPv4Network(user_network, strict=False)
        if options == "1":
            print_network_info(network)
        elif options == "2":
            while True:
                subnets_Needed = input("\nHow many subnets do you need?:\n")
                try:
                    network_based_subnetting(network, subnets_Needed)
                    break
                except ValueError as err:
                    print(err)
                    print("Try again.")
    elif options == "3":
        host_needed = int(input("\nWrite the number of host needed(2-16777216):\n"))
        host_based_subnetting(host_needed)
    elif options == "4":
        break
    else:
        try_again = input("you didn't write a number between 1-3, wanna try again?(y, n):\n").lower()
        if try_again == "n":
            break
