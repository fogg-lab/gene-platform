function generateRandomAlphanumericId(
    length = 8,
    characters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
) {
    return Array.from(
        { length },
        () => characters[Math.floor(Math.random() * characters.length)]
    ).join("");
}

function strReplaceAll(str, search, replace) {
    return str.replace(new RegExp(search, 'g'), replace);
}

export { generateRandomAlphanumericId, strReplaceAll };
