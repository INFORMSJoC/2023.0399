## Instance Format

All instances follow the format below:

- The **first line** contains a single integer **N**, representing the number of items.  
- The **second line** contains a single integer **Q**, representing the total number of distinct classes in the instance.  
- The **third line** contains a single integer **W**, representing the bin capacity.  
- The **fourth line** contains a single integer **C**, representing the maximum number of distinct classes allowed per bin.  
- The next **N** lines each contain two integers:
  - The **weight** of the item  
  - The **class** to which the item belongs (represented by an integer from 0 to Q - 1)

**Note:** Items are **not grouped** by class, size, or any other parameter. Each line corresponds to an individual item with unit demand and an associated class.
