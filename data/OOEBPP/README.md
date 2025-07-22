## Instance Format

All instances follow the **Ordered Open-End Bin Packing Problem (OOEBPP)** format:

- The **first line** contains a single integer **N**, representing the number of items.  
- The **second line** contains two integers:
  - **W** – the bin capacity  
  - A **dummy value** (not used in this problem)

- The next **N** lines each contain three integers:
  - **Item ID** (from 0 to N − 1)  
  - **Weight** of the item  
  - **Priority** value

**Note:**  
Items must be packed in **non-decreasing order of priority**. If multiple items share the same priority, they must respect the order of their **IDs**, with smaller IDs coming first. This rule defines a unique total order over the items.

The priority values were defined this way because the instances were adapted from other problems, and we wanted to preserve the original structure without modifying the data. However, it would be entirely possible to assign a unique priority key to each item and enforce strictly increasing priorities to encode the ordering directly.
